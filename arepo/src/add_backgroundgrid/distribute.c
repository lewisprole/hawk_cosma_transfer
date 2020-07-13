/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/add_backgroundgrid/distribute.c
 * \date        MM/YYYY
 * \author
 * \brief
 * \details
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 */

#include <mpi.h>
#include "../allvars.h"
#include "../proto.h"
#include "add_bggrid.h"

#ifdef ADDBACKGROUNDGRID

static int find_cells_evaluate(int target, int mode, int thread_id);

/* local data structure for collecting particle/cell data that is sent to other processors if needed */
typedef struct
{
  MyDouble Pos[3];
  MyFloat Hsml;
  MyFloat Weight;
  MyFloat Mass;
  MyFloat InternalEnergy;
  MyFloat Momentum[3];
#ifdef MHD
  MyFloat B[3];
#endif
#ifdef METALS
  MyFloat Metallicity;
#endif
#ifdef REFINEMENT_RPS
  MyFloat RPSGalaxyMass;
#endif
  int Firstnode;
} data_in;

static data_in *DataGet;

/*! \brief Routine that fills the relevant particle/cell data into the input
 *         structure defined above. Needed by generic_comm_helpers2.
 *
 *  \param[out] in Data structure to fill.
 *  \param[in] i Index of particle in P and SphP arrays.
 *  \param[in] firstnode First note of communication.
 *
 *  \return void
 */
static void particle2in(data_in *in, int i, int firstnode)
{
  in->Pos[0] = P[i].Pos[0];
  in->Pos[1] = P[i].Pos[1];
  in->Pos[2] = P[i].Pos[2];

  in->Hsml = SphP[i].Hsml;

  in->Weight         = SphP[i].Weight;
  in->Mass           = P[i].Mass;
  in->InternalEnergy = SphP[i].Utherm * P[i].Mass;

  int k;
  for(k = 0; k < 3; k++)
    in->Momentum[k] = P[i].Vel[k] * P[i].Mass;

#ifdef MHD
  for(k = 0; k < 3; k++)
    in->B[k] = SphP[i].B[k];
#endif

#ifdef METALS
  in->Metallicity = SphP[i].Metallicity;
#endif

#ifdef REFINEMENT_RPS
  in->RPSGalaxyMass = SphP[i].RPSGalaxyMass;
#endif

  in->Firstnode = firstnode;
}

/* local data structure that holds results acquired on remote processors */
typedef struct
{
  char nothing;
} data_out;

static data_out *DataResult;

/*! \brief Routine to store or combine result data. Needed by
 *         generic_comm_helpers2.
 *
 *  \param[in] out Data to be moved to appropriate variables in global
 *  particle and cell data arrays (P, SphP,...)
 *  \param[in] i Index of particle in P and SphP arrays
 *  \param[in] mode Mode of function: local particles or information that was
 *  communicated from other tasks and has to be added locally?
 *
 *  \return void
 */
static void out2particle(data_out *out, int i, int mode) { return; }

#include "../generic_comm_helpers2.h"

/*! \brief Routine that defines what to do with local particles.
 *
 *  Calls the *_evaluate function in MODE_LOCAL_PARTICLES.
 *
 *  \return void
 */
static void kernel_local(void)
{
  int idx;
#ifdef GENERIC_ASYNC
  int flag = 0;
#endif
#pragma omp parallel private(i)
  {
    int j, threadid = get_thread_num();
#ifdef GENERIC_ASYNC
    int count = 0;
#endif

    for(j = 0; j < NTask; j++)
      Thread[threadid].Exportflag[j] = -1;

    while(1)
      {
        if(Thread[threadid].ExportSpace < MinSpace)
          break;

#ifdef GENERIC_ASYNC
        if(threadid == 0)
          {
            if((count & POLLINGINTERVAL) == 0)
              if(generic_polling_primary(count, TimeBinsGravity.NActiveParticles))
                flag = 1;

            count++;
          }

        if(flag)
          break;
#endif

#pragma omp atomic capture
        idx = NextParticle++;

        if(idx >= TimeBinsGravity.NActiveParticles)
          break;

        int i = TimeBinsGravity.ActiveParticleList[idx];
        if(i < 0)
          continue;

        find_cells_evaluate(i, MODE_LOCAL_PARTICLES, threadid);
      }
  }
}

/*! \brief Routine that defines what to do with imported particles.
 *
 *  Calls the *_evaluate function in MODE_IMPORTED_PARTICLES.
 *
 *  \return void
 */
static void kernel_imported(void)
{
  /* now do the particles that were sent to us */
  int i, cnt = 0;
#pragma omp parallel private(i)
  {
    int threadid = get_thread_num();
#ifdef GENERIC_ASYNC
    int count = 0;
#endif

    while(1)
      {
#pragma omp atomic capture
        i = cnt++;

        if(i >= Nimport)
          break;

#ifdef GENERIC_ASYNC
        if(threadid == 0)
          {
            if((count & POLLINGINTERVAL) == 0)
              generic_polling_secondary();
          }

        count++;
#endif

        find_cells_evaluate(i, MODE_IMPORTED_PARTICLES, threadid);
      }
  }
}

/*! \brief Main function to distribute hydro quantities over a kernel average.
 *
 *  \return void
 */
void distribute_particles(void)
{
  mpi_printf("ADD BACKGROUND GRID: distributing the fluid quantities\n");

  generic_set_MaxNexport();

  generic_comm_pattern(TimeBinsGravity.NActiveParticles, kernel_local, kernel_imported);

#ifdef MHD
  /* now divide the B field in each cell by the weight (sum of the wk's,
     which we stored in SphP.divB */
  for(int idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      int i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      if(P[i].ID >= IDNew)
        {
          int j;
          if(SphP[i].DivB > 0)
            for(j = 0; j < 3; j++)
              SphP[i].B[j] /= SphP[i].DivB;
        }
    }
#endif

  mpi_printf("ADD BACKGROUND GRID: done\n");
}

int find_cells_evaluate(int target, int mode, int thread_id)
{
  int j, n, numnodes, *firstnode;
  double h, h2, hinv, hinv3;
  MyDouble dx, dy, dz, r;
  MyDouble *pos;
#ifdef PERIODIC
  double xtmp, ytmp, ztmp;
#endif

  data_in local, *target_data;
  data_out out;
  out.nothing = 0;

  if(mode == MODE_LOCAL_PARTICLES)
    {
      particle2in(&local, target, 0);
      target_data = &local;

      numnodes  = 1;
      firstnode = NULL;
    }
  else
    {
      target_data = &DataGet[target];

      generic_get_numnodes(target, &numnodes, &firstnode);
    }

  pos  = target_data->Pos;
  h    = target_data->Hsml;
  h2   = h * h;
  hinv = 1.0 / h;
#ifndef TWODIMS
  hinv3 = hinv * hinv * hinv;
#else
  hinv3 = hinv * hinv / boxSize_Z;
#endif

  int nfound = ngb_treefind_variable_threads(pos, h, target, mode, thread_id, numnodes, firstnode);

  double wsum = 0;

  for(n = 0; n < nfound; n++)
    {
      j = Thread[thread_id].Ngblist[n];

      if(P[j].ID < IDNew)
        continue;

      dx = NGB_PERIODIC_LONG_X(pos[0] - P[j].Pos[0]);
      dy = NGB_PERIODIC_LONG_Y(pos[1] - P[j].Pos[1]);
      dz = NGB_PERIODIC_LONG_Z(pos[2] - P[j].Pos[2]);

      double r2 = dx * dx + dy * dy + dz * dz;

      if(r2 < h2)
        {
          r = sqrt(r2);

          double u = r * hinv;
          double wk;
          if(u < 0.5)
            wk = hinv3 * (KERNEL_COEFF_1 + KERNEL_COEFF_2 * (u - 1) * u * u);
          else
            wk = hinv3 * KERNEL_COEFF_5 * (1.0 - u) * (1.0 - u) * (1.0 - u);

          double weight = SphP[j].Volume * wk / target_data->Weight;

          wsum += weight;

          P[j].Mass += target_data->Mass * weight;
          SphP[j].Energy += target_data->InternalEnergy * weight;

          int k;
          for(k = 0; k < 3; k++)
            SphP[j].Momentum[k] += target_data->Momentum[k] * weight;

#ifdef MHD
          for(k = 0; k < 3; k++)
            SphP[j].B[k] += target_data->B[k] * weight;
          SphP[j].DivB += wk;
#endif
#ifdef METALS
          SphP[j].Metallicity += target_data->Metallicity * target_data->Mass * weight;
#endif
#ifdef REFINEMENT_RPS
          SphP[j].RPSGalaxyMass += target_data->RPSGalaxyMass * weight;
#endif
        }
    }

  if(wsum > 1.01)
    {
      printf("wsum=%g, Weight=%g, target=%d\n", wsum, target_data->Weight, target);
      terminate("bla");
    }

  /* Now collect the result at the right place */
  if(mode == MODE_LOCAL_PARTICLES)
    out2particle(&out, target, MODE_LOCAL_PARTICLES);
  else
    DataResult[target] = out;

  return 0;
}

#endif
