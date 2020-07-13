/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/domain_box.c
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

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>

#include "allvars.h"
#include "domain.h"
#include "proto.h"
#include "voronoi.h"

void domain_displacePosition(double *pos, enum domain_displace_mode mode)
{
  if(mode == DISPLACE_POSITION_FORWARD)
    {
      double xtmp, ytmp, ztmp;
      pos[0] = WRAP_X(pos[0] + All.GlobalDisplacementVector[0]);
      pos[1] = WRAP_Y(pos[1] + All.GlobalDisplacementVector[1]);
      pos[2] = WRAP_Z(pos[2] + All.GlobalDisplacementVector[2]);
    }
  else if(mode == DISPLACE_POSITION_BACKWARD)
    {
      double xtmp, ytmp, ztmp;
      pos[0] = WRAP_X(pos[0] - All.GlobalDisplacementVector[0]);
      pos[1] = WRAP_Y(pos[1] - All.GlobalDisplacementVector[1]);
      pos[2] = WRAP_Z(pos[2] - All.GlobalDisplacementVector[2]);
    }
  else
    terminate("Unkown mode %d.", (int)mode);
}

static void domain_displacePositions(enum domain_displace_mode mode)
{
  for(int i = 0; i < NumPart; i++)
    {
      if(P[i].ID == 0 && P[i].Mass == 0) /* derefined */
        continue;

      domain_displacePosition(P[i].Pos, mode);

      if(i < NumGas)
        domain_displacePosition(SphP[i].Center, mode);

#if defined(BLACK_HOLES) && defined(BH_FRICTION)
      if(P[i].Type == 5)
        {
          domain_displacePosition(BPP(i).BH_MinPotPos_Extended, mode);
          domain_displacePosition(BPP(i).BH_MinPotPos_Previous, mode);
        }
#endif

#if defined(BLACK_HOLES) && defined(MEASURE_POTMIN_AROUND_BH)
      if(P[i].Type == 5)
        {
          domain_displacePosition(BPP(i).BH_MinPotPos, mode);
        }
#endif
    }

#ifdef BH_BASED_CGM_ZOOM
  domain_displacePosition(All.BlackHolePosition, mode);
#endif

#ifdef PLACEHIGHRESREGION
  domain_displacePosition(All.Xmintot[1], mode);
  domain_displacePosition(All.Xmaxtot[1], mode);
  domain_displacePosition(All.Corner[1], mode);
  domain_displacePosition(All.UpperCorner[1], mode);
#endif
}

/*! \brief Finds the extent of the global domain grid.
 *
 *  The minimum extent is the box size.
 *
 *  \return void
 */
void domain_findExtent(void)
{
  int i, j;
  double len, xmin[3], xmax[3], xmin_glob[3], xmax_glob[3];

  /* determine local extension */
  for(j = 0; j < 3; j++)
    {
#ifdef PERIODIC
      /* preset to simulation box */
      xmin[j] = 0;
      xmax[j] = boxSize;
#else
      xmin[j] = MAX_REAL_NUMBER;
      xmax[j] = -MAX_REAL_NUMBER;
#endif
    }
    // Take care of stretched box
#ifdef LONG_X
  xmax[0] = boxSize_X;
#endif
#ifdef LONG_Y
  xmax[1] = boxSize_Y;
#endif
#ifdef LONG_Z
  xmax[2] = boxSize_Z;
#endif

  for(i = 0; i < NumPart; i++)
    {
#ifdef ADDBACKGROUNDGRID
      if(P[i].Type != 0)
        continue;
#endif
      for(j = 0; j < 3; j++)
        {
          if(xmin[j] > P[i].Pos[j])
            xmin[j] = P[i].Pos[j];

          if(xmax[j] < P[i].Pos[j])
            xmax[j] = P[i].Pos[j];
        }
    }

  MPI_Allreduce(xmin, xmin_glob, 3, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(xmax, xmax_glob, 3, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

#ifdef ADDBACKGROUNDGRID
  for(j = 0; j < 3; j++)
    if(xmax_glob[j] < All.BoxSize)
      xmax_glob[j] = All.BoxSize;

  for(j = 0; j < 3; j++)
    if(xmin_glob[j] > 0)
      xmin_glob[j] = 0;
#endif

  len = 0;
  for(j = 0; j < 3; j++)
    if(xmax_glob[j] - xmin_glob[j] > len)
      len = xmax_glob[j] - xmin_glob[j];

#if defined(GRAVITY_NOT_PERIODIC) && !defined(ADDBACKGROUNDGRID)
  len *= 1.2; /* enlarge box a bit to avoid triggering of an out of box recovery */
#else
  len *= 1.00001;
#endif

#ifndef AMR

#if defined(DO_NOT_RANDOMIZE_DOMAINCENTER) || !defined(GRAVITY_NOT_PERIODIC) || defined(ONEDIMS) || defined(TWODIMS)
  for(j = 0; j < 3; j++)
    {
      DomainCenter[j] = 0.5 * (xmin_glob[j] + xmax_glob[j]);
      DomainCorner[j] = 0.5 * (xmin_glob[j] + xmax_glob[j]) - 0.5 * len;
    }
#else
  for(j = 0; j < 3; j++)
    {
      DomainCenter[j] = 0.5 * (xmin_glob[j] + xmax_glob[j]);
      DomainCenter[j] += (2. * get_random_number() - 1.) * 0.5 * len;
    }

  MPI_Bcast(DomainCenter, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  len *= 2;

  for(j = 0; j < 3; j++)
    DomainCorner[j] = DomainCenter[j] - 0.5 * len;
#endif

  DomainLen = len;

#else

  len = fmax(boxSize_X, fmax(boxSize_Y, boxSize_Z));
  DomainLen = len;
  DomainCorner[0] = 0.;
  DomainCorner[1] = 0.;
  DomainCorner[2] = 0.;

  DomainCenter[0] = DomainLen / 2.;
  DomainCenter[1] = DomainLen / 2.;
  DomainCenter[2] = DomainLen / 2.;
#endif

  DomainInverseLen = 1.0 / DomainLen;
  DomainFac        = 1.0 / len * (((peanokey)1) << (BITS_PER_DIMENSION));
  DomainBigFac     = (DomainLen / (((long long)1) << 52));
}

/*! This function makes sure that all particle coordinates (Pos) are
 *  periodically mapped onto the interval [0, BoxSize].  After this function
 *  has been called, a new domain decomposition should be done, which will
 *  also force a new tree construction.
 */
#ifdef PERIODIC
/*! \brief Makes sure all particles are within box.
 *
 *  This function makes sure that all particle coordinates (Pos) are
 *  periodically mapped onto the interval [0, BoxSize].  After this function
 *  has been called, a new domain decomposition should be done, which will
 *  also force a new tree construction.
 *
 *  \return void
 */
void do_box_wrapping(void)
{
  int j;
  double boxsize[3];

#ifdef ADDBACKGROUNDGRID
  return;
#endif

  for(j = 0; j < 3; j++)
    boxsize[j] = All.BoxSize;

#ifdef LONG_X
  boxsize[0] *= LONG_X;
#endif
#ifdef LONG_Y
  boxsize[1] *= LONG_Y;
#endif
#ifdef LONG_Z
  boxsize[2] *= LONG_Z;
#endif

  int i;

#if !defined(GRAVITY_NOT_PERIODIC) && !defined(DO_NOT_RANDOMIZE_DOMAINCENTER) && defined(SELFGRAVITY) && (NUMDIMS > 2)
  domain_displacePositions(DISPLACE_POSITION_BACKWARD);

  if(ThisTask == 0)
    {
      double prefac = 1.;
#ifdef PLACEHIGHRESREGION
      prefac = 0.5;
#endif
      for(j = 0; j < 3; j++)
        All.GlobalDisplacementVector[j] = (get_random_number() - 0.5) * boxsize[j] * prefac;
    }

  mpi_printf("DOMAIN: New global displacement vector: %g, %g, %g, box=%g, %g, %g\n", All.GlobalDisplacementVector[0],
             All.GlobalDisplacementVector[1], All.GlobalDisplacementVector[2], boxsize[0], boxsize[1], boxsize[2]);
  MPI_Bcast(All.GlobalDisplacementVector, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  domain_displacePositions(DISPLACE_POSITION_FORWARD);
#endif

  for(i = 0; i < NumPart; i++)
    {
#if defined(VORONOI_DYNAMIC_UPDATE) || defined(AMR_CONNECTIONS)
      if(i < NumGas)
        trans_table[i].wrapped = 0;
#endif

#if defined(GRAVITY_NOT_PERIODIC)
      if(P[i].Type != 0)
        continue;
#endif

#if !defined(REFLECTIVE_X)
      while(P[i].Pos[0] < 0)
        {
          P[i].Pos[0] += boxsize[0];
#ifdef VORONOI_DYNAMIC_UPDATE
          if(i < NumGas)
            trans_table[i].wrapped |= 1;
#endif
        }

      while(P[i].Pos[0] >= boxsize[0])
        {
          P[i].Pos[0] -= boxsize[0];
#ifdef VORONOI_DYNAMIC_UPDATE
          if(i < NumGas)
            trans_table[i].wrapped |= 2;
#endif
        }
#else
      if(P[i].Pos[0] < 0 || P[i].Pos[0] >= boxsize[0])
        terminate("i=%d ID=%" MYIDTYPE_PRI " type=%d moved out of box. x=%g", i, P[i].ID, P[i].Type, P[i].Pos[0]);
#endif

#if !defined(REFLECTIVE_Y)
      while(P[i].Pos[1] < 0)
        {
          P[i].Pos[1] += boxsize[1];
#ifdef VORONOI_DYNAMIC_UPDATE
          if(i < NumGas)
            trans_table[i].wrapped |= 4;
#endif
        }

      while(P[i].Pos[1] >= boxsize[1])
        {
          P[i].Pos[1] -= boxsize[1];
#ifdef VORONOI_DYNAMIC_UPDATE
          if(i < NumGas)
            trans_table[i].wrapped |= 8;
#endif
        }

#else
      if(P[i].Pos[1] < 0 || P[i].Pos[1] >= boxsize[1])
        terminate("i=%d ID=%" MYIDTYPE_PRI " type=%d moved out of box. y=%g", i, P[i].ID, P[i].Type, P[i].Pos[1]);
#endif

#if !defined(REFLECTIVE_Z)
      while(P[i].Pos[2] < 0)
        {
          P[i].Pos[2] += boxsize[2];
#ifdef VORONOI_DYNAMIC_UPDATE
          if(i < NumGas)
            trans_table[i].wrapped |= 16;
#endif
        }

      while(P[i].Pos[2] >= boxsize[2])
        {
          P[i].Pos[2] -= boxsize[2];
#ifdef VORONOI_DYNAMIC_UPDATE
          if(i < NumGas)
            trans_table[i].wrapped |= 32;
#endif
        }
#else
      if(P[i].Pos[2] < 0 || P[i].Pos[2] >= boxsize[2])
        terminate("i=%d ID=%" MYIDTYPE_PRI " type=%d moved out of box. z=%g", i, P[i].ID, P[i].Type, P[i].Pos[2]);
#endif
    }
}
#endif
