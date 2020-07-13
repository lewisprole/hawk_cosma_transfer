/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/parameters.c
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

#include <gsl/gsl_rng.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include "allvars.h"
#include "proto.h"

/*! \file parameters.c
 *  \brief  Parses the parameter file
 *
 *  This file contains the routine to parse the parameter file.
 *  Additionally the output list is also parsed.
 */

/*! \brief This function parses the parameter file.
 *
 *  Each parameter is defined by a keyword (`tag'), and can be either
 *  of type douple, int, or character string. Three arrays containing the name,
 *  type and address of the parameter are filled first. The routine then parses
 *  the parameter file and fills the referenced variables. The routine makes
 *  sure that each parameter appears exactly once in the parameter file,
 *  otherwise error messages are produced that complain about the missing
 *  parameters.
 *
 *  \param[in] fname The file name of the parameter file
 *
 *  \return void
 */
void read_parameter_file(char *fname)
{
#define REAL 1
#define STRING 2
#define INT 3

  FILE *fd, *fdout;
  char buf[MAXLEN_PARAM_TAG + MAXLEN_PARAM_VALUE + 200], buf1[MAXLEN_PARAM_TAG + 200], buf2[MAXLEN_PARAM_VALUE + 200],
      buf3[MAXLEN_PARAM_TAG + MAXLEN_PARAM_VALUE + 400];
  int i, j, nt;
  int id[MAX_PARAMETERS];
  void *addr[MAX_PARAMETERS];
  char tag[MAX_PARAMETERS][MAXLEN_PARAM_TAG];
  int param_handled[MAX_PARAMETERS];
  int errorFlag = 0;

  All.StarformationOn = 0; /* defaults */

  for(i = 0; i < MAX_PARAMETERS; i++)
    {
      param_handled[i] = 0;
    }

  if(sizeof(long long) != 8)
    {
      mpi_terminate("\nType `long long' is not 64 bit on this platform. Stopping.\n\n");
    }

  if(sizeof(int) != 4)
    {
      mpi_terminate("\nType `int' is not 32 bit on this platform. Stopping.\n\n");
    }

  if(sizeof(float) != 4)
    {
      mpi_terminate("\nType `float' is not 32 bit on this platform. Stopping.\n\n");
    }

  if(sizeof(double) != 8)
    {
      mpi_terminate("\nType `double' is not 64 bit on this platform. Stopping.\n\n");
    }

  if(ThisTask == 0) /* read parameter file on process 0 */
    {
      nt = 0;

      strcpy(tag[nt], "InitCondFile");
      addr[nt] = All.InitCondFile;
      id[nt++] = STRING;

      strcpy(tag[nt], "OutputDir");
      addr[nt] = All.OutputDir;
      id[nt++] = STRING;

#ifdef TOLERATE_WRITE_ERROR
      strcpy(tag[nt], "AlternativeOutputDir");
      addr[nt] = AlternativeOutputDir;
      id[nt++] = STRING;
#endif

      strcpy(tag[nt], "SnapshotFileBase");
      addr[nt] = All.SnapshotFileBase;
      id[nt++] = STRING;

      strcpy(tag[nt], "ResubmitCommand");
      addr[nt] = All.ResubmitCommand;
      id[nt++] = STRING;

      strcpy(tag[nt], "OutputListFilename");
      addr[nt] = All.OutputListFilename;
      id[nt++] = STRING;

      strcpy(tag[nt], "OutputListOn");
      addr[nt] = &All.OutputListOn;
      id[nt++] = INT;

      strcpy(tag[nt], "Omega0");
      addr[nt] = &All.Omega0;
      id[nt++] = REAL;

      strcpy(tag[nt], "OmegaBaryon");
      addr[nt] = &All.OmegaBaryon;
      id[nt++] = REAL;

      strcpy(tag[nt], "OmegaLambda");
      addr[nt] = &All.OmegaLambda;
      id[nt++] = REAL;

      strcpy(tag[nt], "HubbleParam");
      addr[nt] = &All.HubbleParam;
      id[nt++] = REAL;

      strcpy(tag[nt], "BoxSize");
      addr[nt] = &All.BoxSize;
      id[nt++] = REAL;

      strcpy(tag[nt], "PeriodicBoundariesOn");
      addr[nt] = &All.PeriodicBoundariesOn;
      id[nt++] = INT;

      strcpy(tag[nt], "MaxMemSize");
      addr[nt] = &All.MaxMemSize;
      id[nt++] = INT;

      strcpy(tag[nt], "TimeOfFirstSnapshot");
      addr[nt] = &All.TimeOfFirstSnapshot;
      id[nt++] = REAL;

      strcpy(tag[nt], "CpuTimeBetRestartFile");
      addr[nt] = &All.CpuTimeBetRestartFile;
      id[nt++] = REAL;

#ifdef REDUCE_FLUSH
      strcpy(tag[nt], "FlushCpuTimeDiff");
      addr[nt] = &All.FlushCpuTimeDiff;
      id[nt++] = REAL;
#endif

      strcpy(tag[nt], "TimeBetStatistics");
      addr[nt] = &All.TimeBetStatistics;
      id[nt++] = REAL;

      strcpy(tag[nt], "TimeBegin");
      addr[nt] = &All.TimeBegin;
      id[nt++] = REAL;

      strcpy(tag[nt], "TimeMax");
      addr[nt] = &All.TimeMax;
      id[nt++] = REAL;

      strcpy(tag[nt], "TimeBetSnapshot");
      addr[nt] = &All.TimeBetSnapshot;
      id[nt++] = REAL;

      strcpy(tag[nt], "UnitVelocity_in_cm_per_s");
      addr[nt] = &All.UnitVelocity_in_cm_per_s;
      id[nt++] = REAL;

      strcpy(tag[nt], "UnitLength_in_cm");
      addr[nt] = &All.UnitLength_in_cm;
      id[nt++] = REAL;

      strcpy(tag[nt], "UnitMass_in_g");
      addr[nt] = &All.UnitMass_in_g;
      id[nt++] = REAL;

      strcpy(tag[nt], "ErrTolIntAccuracy");
      addr[nt] = &All.ErrTolIntAccuracy;
      id[nt++] = REAL;

      strcpy(tag[nt], "ErrTolTheta");
      addr[nt] = &All.ErrTolTheta;
      id[nt++] = REAL;

      strcpy(tag[nt], "ErrTolForceAcc");
      addr[nt] = &All.ErrTolForceAcc;
      id[nt++] = REAL;

      strcpy(tag[nt], "MaxSizeTimestep");
      addr[nt] = &All.MaxSizeTimestep;
      id[nt++] = REAL;

      strcpy(tag[nt], "MinSizeTimestep");
      addr[nt] = &All.MinSizeTimestep;
      id[nt++] = REAL;

#ifdef LEGACY_DISPLACEMENT_CONSTRAINT
      strcpy(tag[nt], "MaxRMSDisplacementFac");
      addr[nt] = &All.MaxRMSDisplacementFac;
      id[nt++] = REAL;
#endif

      strcpy(tag[nt], "CourantFac");
      addr[nt] = &All.CourantFac;
      id[nt++] = REAL;

      strcpy(tag[nt], "LimitUBelowThisDensity");
      addr[nt] = &All.LimitUBelowThisDensity;
      id[nt++] = REAL;

      strcpy(tag[nt], "LimitUBelowCertainDensityToThisValue");
      addr[nt] = &All.LimitUBelowCertainDensityToThisValue;
      id[nt++] = REAL;

      strcpy(tag[nt], "DesNumNgb");
      addr[nt] = &All.DesNumNgb;
      id[nt++] = INT;

      strcpy(tag[nt], "MultipleDomains");
      addr[nt] = &All.MultipleDomains;
      id[nt++] = INT;

      strcpy(tag[nt], "TopNodeFactor");
      addr[nt] = &All.TopNodeFactor;
      id[nt++] = REAL;

      strcpy(tag[nt], "ActivePartFracForNewDomainDecomp");
      addr[nt] = &All.ActivePartFracForNewDomainDecomp;
      id[nt++] = REAL;

#ifdef BECDM
      strcpy(tag[nt], "AxionMassEv");
      addr[nt] = &All.AxionMassEv;
      id[nt++] = REAL;
#endif

#ifdef SUBFIND
      strcpy(tag[nt], "DesLinkNgb");
      addr[nt] = &All.DesLinkNgb;
      id[nt++] = INT;

      strcpy(tag[nt], "ErrTolThetaSubfind");
      addr[nt] = &All.ErrTolThetaSubfind;
      id[nt++] = REAL;
#endif

#ifdef WINDTUNNEL
      strcpy(tag[nt], "InjectionDensity");
      addr[nt] = &All.InjectionDensity;
      id[nt++] = REAL;

      strcpy(tag[nt], "InjectionVelocity");
      addr[nt] = &All.InjectionVelocity;
      id[nt++] = REAL;

      strcpy(tag[nt], "InjectionUtherm");
      addr[nt] = &All.InjectionUtherm;
      id[nt++] = REAL;

      strcpy(tag[nt], "InjectionRegion");
      addr[nt] = &All.InjectionRegion;
      id[nt++] = REAL;

      strcpy(tag[nt], "InjectionVolume");
      addr[nt] = &All.InjectionVolume;
      id[nt++] = REAL;

#ifdef WINDTUNNEL_EXTERNAL_SOURCE
      strcpy(tag[nt], "WindTunnelExternalSourceFile");
      addr[nt] = All.WindTunnelExternalSourceFile;
      id[nt++] = STRING;
#endif

#if defined(WINDTUNNEL_FIXVARIABLESININJECTIONREGION) && defined(MHD)
      strcpy(tag[nt], "InjectionBx_InGauss");
      addr[nt] = &All.InjectionBx_InGauss;
      id[nt++] = REAL;

      strcpy(tag[nt], "InjectionBy_InGauss");
      addr[nt] = &All.InjectionBy_InGauss;
      id[nt++] = REAL;

      strcpy(tag[nt], "InjectionBz_InGauss");
      addr[nt] = &All.InjectionBz_InGauss;
      id[nt++] = REAL;
#endif

#if defined(WINDTUNNEL) && defined(WINDTUNNEL_READ_IN_BFIELD)
      strcpy(tag[nt], "WindtunnelReadIn_InputFileName");
      addr[nt] = All.WindtunnelReadIn_InputFileName;
      id[nt++] = STRING;
#endif

#endif /* WINDTUNNEL */

#ifdef CIRCUMSTELLAR
#ifdef CIRCUMSTELLAR_REFINEMENTS
      strcpy(tag[nt], "CircumstellarDerefinementDistance");
      addr[nt] = &All.CircumstellarDerefinementDistance;
      id[nt++] = REAL;
#endif
#ifdef CIRCUMSTELLAR_IRRADIATION
      strcpy(tag[nt], "IrradiationTempScaling");
      addr[nt] = &All.IrradiationTempScaling;
      id[nt++] = REAL;
#endif
#endif

#if defined(LOCALLY_ISOTHERM_DISK)
      strcpy(tag[nt], "AspectRatio");
      addr[nt] = &All.AspectRatio;
      id[nt++] = REAL;

      strcpy(tag[nt], "CentralMass");
      addr[nt] = &All.CentralMass;
      id[nt++] = REAL;

      strcpy(tag[nt], "InnerRadius");
      addr[nt] = &All.InnerRadius;
      id[nt++] = REAL;

      strcpy(tag[nt], "OuterRadius");
      addr[nt] = &All.OuterRadius;
      id[nt++] = REAL;
#endif /* #if defined(LOCALLY_ISOTHERM_DISK) */

#ifdef CIRCUMSTELLAR_WBOUNDARIES
      strcpy(tag[nt], "InnerRadius");
      addr[nt] = &All.inner_radius;
      id[nt++] = REAL;

      strcpy(tag[nt], "OuterRadius");
      addr[nt] = &All.outer_radius;
      id[nt++] = REAL;

      strcpy(tag[nt], "CircumstellarBoundaryDensity");
      addr[nt] = &All.CircumstellarBoundaryDensity;
      id[nt++] = REAL;

      strcpy(tag[nt], "EvanescentBoundaryStrength");
      addr[nt] = &All.EvanescentBoundaryStrength;
      id[nt++] = REAL;
#endif

#ifdef CENTRAL_MASS_POTENTIAL
      strcpy(tag[nt], "CentralMass");
      addr[nt] = &All.CentralMass;
      id[nt++] = REAL;

      strcpy(tag[nt], "SofteningCentral");
      addr[nt] = &All.SofteningCentral;
      id[nt++] = REAL;
#endif

#if defined(ACCRETE_ONTO_CENTRAL_POTENTIAL) && defined(CENTRAL_MASS_POTENTIAL)
      strcpy(tag[nt], "CentralAccretionRadius");
      addr[nt] = &All.CentralAccretionRadius;
      id[nt++] = REAL;
#endif

#ifdef PERTURB_VELOCITIES
      strcpy(tag[nt], "VelocityPerturbation");
      addr[nt] = &All.VelocityPerturbation;
      id[nt++] = REAL;
#endif

#ifdef STAR_PLANET_POTENTIAL
      strcpy(tag[nt], "MassRatio");
      addr[nt] = &All.MassRatio;
      id[nt++] = REAL;

      strcpy(tag[nt], "SofteningPlanet");
      addr[nt] = &All.SofteningPlanet;
      id[nt++] = REAL;

      strcpy(tag[nt], "PlanetGrowthTime");
      addr[nt] = &All.PlanetGrowthTime;
      id[nt++] = REAL;
#endif

#ifdef BINARY_POTENTIAL
      strcpy(tag[nt], "BinaryMassRatio");
      addr[nt] = &All.BinaryMassRatio;
      id[nt++] = REAL;

      strcpy(tag[nt], "BinarySoftening");
      addr[nt] = &All.BinarySoftening;
      id[nt++] = REAL;

      strcpy(tag[nt], "BinaryGrowthTime");
      addr[nt] = &All.BinaryGrowthTime;
      id[nt++] = REAL;

      strcpy(tag[nt], "BinaryEccentricity");
      addr[nt] = &All.BinaryEccentricity;
      id[nt++] = REAL;

      strcpy(tag[nt], "BinaryBarycentricCoord");
      addr[nt] = &All.BinaryBarycentricCoord;
      id[nt++] = INT;
#endif

#if defined(ISOTHERM_EQS) || defined(VS_TURB) || defined(AB_TURB)
      strcpy(tag[nt], "IsoSoundSpeed");
      addr[nt] = &All.IsoSoundSpeed;
      id[nt++] = REAL;
#endif

#if(defined(VS_TURB) || defined(AB_TURB)) && defined(POWERSPEC_GRID)
      strcpy(tag[nt], "TimeBetTurbSpectrum");
      addr[nt] = &All.TimeBetTurbSpectrum;
      id[nt++] = REAL;
#endif

#if defined(AB_TURB) && defined(AB_TURB_DECAYING)
      strcpy(tag[nt], "TimeTurbDecay");
      addr[nt] = &All.TimeTurbDecay;
      id[nt++] = REAL;
#endif

#ifdef DVR_RENDER
#if(DVR_RENDER == 1)
      strcpy(tag[nt], "DvrRenderTimeIntverallInGyr");
      addr[nt] = &All.DvrRenderTimeIntverallInGyr;
      id[nt++] = REAL;

      strcpy(tag[nt], "DvrPixelX");
      addr[nt] = &All.DvrPixelX;
      id[nt++] = INT;

      strcpy(tag[nt], "DvrPixelY");
      addr[nt] = &All.DvrPixelY;
      id[nt++] = INT;
#endif
      strcpy(tag[nt], "DvrTauScaleFactor");
      addr[nt] = &All.DvrTauScaleFactor;
      id[nt++] = REAL;

      strcpy(tag[nt], "DvrTauFloor");
      addr[nt] = &All.DvrTauFloor;
      id[nt++] = REAL;

      strcpy(tag[nt], "DvrOutputDir");
      addr[nt] = &All.DvrOutputDir;
      id[nt++] = STRING;
#endif
      strcpy(tag[nt], "MaxNumNgbDeviation");
      addr[nt] = &All.MaxNumNgbDeviation;
      id[nt++] = REAL;

      strcpy(tag[nt], "ComovingIntegrationOn");
      addr[nt] = &All.ComovingIntegrationOn;
      id[nt++] = INT;

      strcpy(tag[nt], "ICFormat");
      addr[nt] = &All.ICFormat;
      id[nt++] = INT;

      strcpy(tag[nt], "SnapFormat");
      addr[nt] = &All.SnapFormat;
      id[nt++] = INT;

      strcpy(tag[nt], "NumFilesPerSnapshot");
      addr[nt] = &All.NumFilesPerSnapshot;
      id[nt++] = INT;

      strcpy(tag[nt], "NumFilesWrittenInParallel");
      addr[nt] = &All.NumFilesWrittenInParallel;
      id[nt++] = INT;

      strcpy(tag[nt], "ResubmitOn");
      addr[nt] = &All.ResubmitOn;
      id[nt++] = INT;

      strcpy(tag[nt], "CoolingOn");
      addr[nt] = &All.CoolingOn;
      id[nt++] = INT;

      strcpy(tag[nt], "StarformationOn");
      addr[nt] = &All.StarformationOn;
      id[nt++] = INT;

      strcpy(tag[nt], "TypeOfTimestepCriterion");
      addr[nt] = &All.TypeOfTimestepCriterion;
      id[nt++] = INT;

      strcpy(tag[nt], "TypeOfOpeningCriterion");
      addr[nt] = &All.TypeOfOpeningCriterion;
      id[nt++] = INT;

      strcpy(tag[nt], "TimeLimitCPU");
      addr[nt] = &All.TimeLimitCPU;
      id[nt++] = REAL;

      strcpy(tag[nt], "GasSoftFactor");
      addr[nt] = &All.GasSoftFactor;
      id[nt++] = REAL;

      for(i = 0; i < NSOFTTYPES; i++)
        {
          char buf[100];
          sprintf(buf, "SofteningComovingType%d", i);
          strcpy(tag[nt], buf);
          addr[nt] = &All.SofteningComoving[i];
          id[nt++] = REAL;
        }

      for(i = 0; i < NSOFTTYPES; i++)
        {
          char buf[100];
          sprintf(buf, "SofteningMaxPhysType%d", i);
          strcpy(tag[nt], buf);
          addr[nt] = &All.SofteningMaxPhys[i];
          id[nt++] = REAL;
        }

      for(i = 0; i < NTYPES; i++)
        {
          char buf[100];
          sprintf(buf, "SofteningTypeOfPartType%d", i);
          strcpy(tag[nt], buf);
          addr[nt] = &All.SofteningTypeOfPartType[i];
          id[nt++] = INT;
        }

#ifdef ADAPTIVE_HYDRO_SOFTENING
      strcpy(tag[nt], "MinimumComovingHydroSoftening");
      addr[nt] = &All.MinimumComovingHydroSoftening;
      id[nt++] = REAL;

      strcpy(tag[nt], "AdaptiveHydroSofteningSpacing");
      addr[nt] = &All.AdaptiveHydroSofteningSpacing;
      id[nt++] = REAL;
#endif

      strcpy(tag[nt], "GravityConstantInternal");
      addr[nt] = &All.GravityConstantInternal;
      id[nt++] = REAL;

      strcpy(tag[nt], "InitGasTemp");
      addr[nt] = &All.InitGasTemp;
      id[nt++] = REAL;

      strcpy(tag[nt], "MinGasTemp");
      addr[nt] = &All.MinGasTemp;
      id[nt++] = REAL;

      strcpy(tag[nt], "MinEgySpec");
      addr[nt] = &All.MinEgySpec;
      id[nt++] = REAL;

      strcpy(tag[nt], "MinimumDensityOnStartUp");
      addr[nt] = &All.MinimumDensityOnStartUp;
      id[nt++] = REAL;

#ifdef NODEREFINE_BACKGROUND_GRID
      strcpy(tag[nt], "MeanVolume");
      addr[nt] = &All.MeanVolume;
      id[nt++] = REAL;
#endif

#ifdef REFINEMENT_CGM
      strcpy(tag[nt], "TargetVolumeRelativeToSFRThreshold");
      addr[nt] = &All.TargetVolumeRelativeToSFRThreshold;
      id[nt++] = REAL;

      strcpy(tag[nt], "MinMassForCGMRefinement");
      addr[nt] = &All.MinMassForCGMRefinement;
      id[nt++] = REAL;

      strcpy(tag[nt], "FracRadiusForCGMRefinement");
      addr[nt] = &All.FracRadiusForCGMRefinement;
      id[nt++] = REAL;
#endif

#ifdef BH_NEW_CENTERING
      strcpy(tag[nt], "BlackHoleCenteringMassMultiplier");
      addr[nt] = &All.BlackHoleCenteringMassMultiplier;
      id[nt++] = REAL;
#endif

#if defined(BLACK_HOLES) && defined(BH_FRICTION)
      strcpy(tag[nt], "BHFrictionCoefficient");
      addr[nt] = &All.BHFrictionCoefficient;
      id[nt++] = REAL;

      strcpy(tag[nt], "BHFrictionAvgTime");
      addr[nt] = &All.BHFrictionAvgTime;
      id[nt++] = REAL;
#endif

#if defined(BLACK_HOLES) && defined(BH_SPIN_EVOLUTION)
      strcpy(tag[nt], "BHInitialSpin");
      addr[nt] = &All.BHInitialSpin;
      id[nt++] = REAL;

      strcpy(tag[nt], "ShakuraSunyaevParameter");
      addr[nt] = &All.ShakuraSunyaevParameter;
      id[nt++] = REAL;
#endif

#ifdef REFINEMENT_AROUND_BH
#ifdef REFINEMENT_AROUND_BH_FIXED
      strcpy(tag[nt], "RefBHRadius");
      addr[nt] = &All.RefBHRadius;
      id[nt++] = REAL;

      strcpy(tag[nt], "RefBHMinCellRadius");
      addr[nt] = &All.RefBHMinCellRadius;
      id[nt++] = REAL;

      strcpy(tag[nt], "RefBHMaxCellRadius");
      addr[nt] = &All.RefBHMaxCellRadius;
      id[nt++] = REAL;
#else
      strcpy(tag[nt], "RefBHRadiusHSML");
      addr[nt] = &All.RefBHRadiusHSML;
      id[nt++] = REAL;

      strcpy(tag[nt], "RefBHMinCellRadiusRBondi");
      addr[nt] = &All.RefBHMinCellRadiusRBondi;
      id[nt++] = REAL;

      strcpy(tag[nt], "RefBHMaxCellRadiusHSML");
      addr[nt] = &All.RefBHMaxCellRadiusHSML;
      id[nt++] = REAL;
#endif
      strcpy(tag[nt], "RefBHMinCellMass");
      addr[nt] = &All.RefBHMinCellMass;
      id[nt++] = REAL;

      strcpy(tag[nt], "RefBHLowerFactorC");
      addr[nt] = &All.RefBHLowerFactorC;
      id[nt++] = REAL;
#endif

#ifdef REFINEMENT_AROUND_DM
      strcpy(tag[nt], "RefinementCellsPerSoftening");
      addr[nt] = &All.RefinementCellsPerSoftening;
      id[nt++] = REAL;
#endif

#ifdef BH_BIPOLAR_FEEDBACK
      strcpy(tag[nt], "BHBipolarTheta");
      addr[nt] = &All.BHBipolarTheta;
      id[nt++] = REAL;

      strcpy(tag[nt], "BHBipolarEfficiency");
      addr[nt] = &All.BHBipolarEfficiency;
      id[nt++] = REAL;

      strcpy(tag[nt], "BHBipolarColdFraction");
      addr[nt] = &All.BHBipolarColdFraction;
      id[nt++] = REAL;

      strcpy(tag[nt], "BHBipolarColdTemp");
      addr[nt] = &All.BHBipolarColdTemp;
      id[nt++] = REAL;
#endif

#ifdef MIN_METALLICITY_ON_STARTUP
      strcpy(tag[nt], "MinimumMetallicityOnStartUp");
      addr[nt] = &All.MinimumMetallicityOnStartUp;
      id[nt++] = REAL;
#endif

#ifndef VORONOI_STATIC_MESH
#ifdef REGULARIZE_MESH_FACE_ANGLE
      strcpy(tag[nt], "CellMaxAngleFactor");
      addr[nt] = &All.CellMaxAngleFactor;
      id[nt++] = REAL;
#else
      strcpy(tag[nt], "CellShapingFactor");
      addr[nt] = &All.CellShapingFactor;
      id[nt++] = REAL;
#endif

      strcpy(tag[nt], "CellShapingSpeed");
      addr[nt] = &All.CellShapingSpeed;
      id[nt++] = REAL;
#endif

#if defined(VORONOI_IMAGES_FOREACHSNAPSHOT) || defined(VORONOI_FREQUENT_IMAGES)
      strcpy(tag[nt], "PicXpixels");
      addr[nt] = &All.PicXpixels;
      id[nt++] = INT;

      strcpy(tag[nt], "PicYpixels");
      addr[nt] = &All.PicYpixels;
      id[nt++] = INT;

      strcpy(tag[nt], "PicXaxis");
      addr[nt] = &All.PicXaxis;
      id[nt++] = INT;

      strcpy(tag[nt], "PicYaxis");
      addr[nt] = &All.PicYaxis;
      id[nt++] = INT;

      strcpy(tag[nt], "PicZaxis");
      addr[nt] = &All.PicZaxis;
      id[nt++] = INT;

      strcpy(tag[nt], "PicXmin");
      addr[nt] = &All.PicXmin;
      id[nt++] = REAL;

      strcpy(tag[nt], "PicXmax");
      addr[nt] = &All.PicXmax;
      id[nt++] = REAL;

      strcpy(tag[nt], "PicYmin");
      addr[nt] = &All.PicYmin;
      id[nt++] = REAL;

      strcpy(tag[nt], "PicYmax");
      addr[nt] = &All.PicYmax;
      id[nt++] = REAL;

      strcpy(tag[nt], "PicZmin");
      addr[nt] = &All.PicZmin;
      id[nt++] = REAL;

      strcpy(tag[nt], "PicZmax");
      addr[nt] = &All.PicZmax;
      id[nt++] = REAL;
#endif

#ifdef VORONOI_FREQUENT_IMAGES
      strcpy(tag[nt], "TimeBetweenImages");
      addr[nt] = &All.TimeBetweenImages;
      id[nt++] = REAL;
#endif

#ifdef SUBBOX_SNAPSHOTS
      strcpy(tag[nt], "SubboxCoordinatesPath");
      addr[nt] = &All.SubboxCoordinatesPath;
      id[nt++] = STRING;

      strcpy(tag[nt], "SubboxMinTime");
      addr[nt] = &All.SubboxMinTime;
      id[nt++] = REAL;

      strcpy(tag[nt], "SubboxMaxTime");
      addr[nt] = &All.SubboxMaxTime;
      id[nt++] = REAL;

      strcpy(tag[nt], "SubboxSyncModulo");
      addr[nt] = &All.SubboxSyncModulo;
      id[nt++] = INT;

      strcpy(tag[nt], "SubboxNumFilesPerSnapshot");
      addr[nt] = &All.SubboxNumFilesPerSnapshot;
      id[nt++] = INT;

      strcpy(tag[nt], "SubboxNumFilesWrittenInParallel");
      addr[nt] = &All.SubboxNumFilesWrittenInParallel;
      id[nt++] = INT;
#endif

#if defined(CIRCUMSTELLAR) && defined(CIRCUMSTELLAR_SINKS)
      strcpy(tag[nt], "CircumstellarSinkRadius");
      addr[nt] = &All.CircumstellarSinkRadius;
      id[nt++] = REAL;
#endif

#if defined(FOF) && (defined(BLACK_HOLES) || defined(GFM_WINDS_VARIABLE) || defined(GFM_BIPOLAR_WINDS) || defined(GFM_WINDS_LOCAL))
      strcpy(tag[nt], "TimeBetOnTheFlyFoF");
      addr[nt] = &All.TimeBetOnTheFlyFoF;
      id[nt++] = REAL;
#endif

#ifdef BLACK_HOLES

#ifdef BH_RELATIVE_NGB_DEVIATION
      strcpy(tag[nt], "DesNumNgbBlackHoleRelDeviationFactor");
      addr[nt] = &All.DesNumNgbBlackHoleRelDeviationFactor;
      id[nt++] = REAL;
#endif

#ifdef BH_THERMALFEEDBACK_ACC
      strcpy(tag[nt], "BlackholeDeltaTemp");
      addr[nt] = &All.BlackholeDeltaTemp;
      id[nt++] = REAL;

      strcpy(tag[nt], "BlackholeDeltaTime");
      addr[nt] = &All.BlackholeDeltaTime;
      id[nt++] = REAL;
#endif

#ifdef BH_BONDI_DENSITY
      strcpy(tag[nt], "BlackHoleAccretionSlope");
      addr[nt] = &All.BlackHoleAccretionSlope;
      id[nt++] = REAL;
#else
      strcpy(tag[nt], "BlackHoleAccretionFactor");
      addr[nt]           = &All.BlackHoleAccretionFactor;
      id[nt++]           = REAL;
#endif

#ifdef BH_BONDI_DISK_VORTICITY
      strcpy(tag[nt], "DiskVorticityRadius");
      addr[nt] = &All.DiskVorticityRadius;
      id[nt++] = REAL;
#endif

      strcpy(tag[nt], "BlackHoleFeedbackFactor");
      addr[nt] = &All.BlackHoleFeedbackFactor;
      id[nt++] = REAL;

      strcpy(tag[nt], "BlackHoleEddingtonFactor");
      addr[nt] = &All.BlackHoleEddingtonFactor;
      id[nt++] = REAL;

      strcpy(tag[nt], "SeedBlackHoleMass");
      addr[nt] = &All.SeedBlackHoleMass;
      id[nt++] = REAL;
#ifdef FOF
      strcpy(tag[nt], "MinFoFMassForNewSeed");
      addr[nt] = &All.MinFoFMassForNewSeed;
      id[nt++] = REAL;
#endif
#ifdef MASSIVE_SEEDS
      strcpy(tag[nt], "DesNumNgbSeed");
      addr[nt] = &All.DesNumNgbSeed;
      id[nt++] = REAL;

      strcpy(tag[nt], "SeedMaxAccretionRadius");
      addr[nt] = &All.SeedMaxAccretionRadius;
      id[nt++] = REAL;
#endif

      strcpy(tag[nt], "DesNumNgbBlackHole");
      addr[nt] = &All.DesNumNgbBlackHole;
      id[nt++] = INT;

      strcpy(tag[nt], "BlackHoleMaxAccretionRadius");
      addr[nt] = &All.BlackHoleMaxAccretionRadius;
      id[nt++] = REAL;

      strcpy(tag[nt], "BlackHoleRadiativeEfficiency");
      addr[nt] = &All.BlackHoleRadiativeEfficiency;
      id[nt++] = REAL;

#ifdef BH_NF_RADIO
      strcpy(tag[nt], "RadioModeMachnumber");
      addr[nt] = &All.RadioModeMachnumber;
      id[nt++] = REAL;

      strcpy(tag[nt], "RadioRelativeBubbleSize");
      addr[nt] = &All.RadioRelativeBubbleSize;
      id[nt++] = REAL;

      strcpy(tag[nt], "RadioRelativeBubbleEnergy");
      addr[nt] = &All.RadioRelativeBubbleEnergy;
      id[nt++] = REAL;

      strcpy(tag[nt], "RadioRelativeMaxDist");
      addr[nt] = &All.RadioRelativeMaxDist;
      id[nt++] = REAL;

      strcpy(tag[nt], "RadioModeMetallicityInSolar");
      addr[nt] = &All.RadioModeMetallicityInSolar;
      id[nt++] = REAL;
#endif

#ifdef BH_BASED_CGM_ZOOM
      strcpy(tag[nt], "CGM_MinRadius");
      addr[nt] = &All.CGM_MinRadius;
      id[nt++] = REAL;

      strcpy(tag[nt], "CGM_MaxRadius");
      addr[nt] = &All.CGM_MaxRadius;
      id[nt++] = REAL;

      strcpy(tag[nt], "IGM_Radius");
      addr[nt] = &All.IGM_Radius;
      id[nt++] = REAL;

      strcpy(tag[nt], "CGM_RefinementFactor");
      addr[nt] = &All.CGM_RefinementFactor;
      id[nt++] = REAL;

      strcpy(tag[nt], "CGM_RadiiRedshift");
      addr[nt] = &All.CGM_RadiiRedshift;
      id[nt++] = REAL;
#endif

#endif /* BLACK_HOLES */

#ifdef MRT

#ifdef MRT_SUBCYCLE
      strcpy(tag[nt], "RTNumSubCycles");
      addr[nt] = &All.RTNumSubCycles;
      id[nt++] = INT;
#else
      All.RTNumSubCycles = 1.0;
#endif

#ifdef MRT_CHEM_SG
      strcpy(tag[nt], "UVKappa");
      addr[nt] = &All.UVKappa;
      id[nt++] = REAL;
#ifdef MRT_IR
      strcpy(tag[nt], "IRKappa");
      addr[nt] = &All.IRKappa;
      id[nt++] = REAL;
#endif
#endif

#ifdef MRT_RIEMANN_HLLE
      strcpy(tag[nt], "HLLEFile");
      addr[nt] = &All.HLLEFile;
      id[nt++] = STRING;
#endif

#ifdef MRT_SINGLE_STAR
      strcpy(tag[nt], "StellarParamFile");
      addr[nt] = &All.StellarParamFile;
      id[nt++] = STRING;
#endif

#ifdef MRT_IR_GRAIN_KAPPA
      strcpy(tag[nt], "GrainKappaPath");
      addr[nt] = &All.GrainKappaPath;
      id[nt++] = STRING;
#endif
#ifdef MRT_BH
#ifdef MRT_BH_PULSED
      strcpy(tag[nt], "AGNPulseTime");
      addr[nt] = &All.AGNPulseTime;
      id[nt++] = REAL;
#endif
#ifdef MRT_BH_BIPOLAR
      strcpy(tag[nt], "MRTBH_OpeningAngle");
      addr[nt] = &All.MRTBH_OpeningAngle;
      id[nt++] = REAL;
#endif

      strcpy(tag[nt], "LogAGNLuminosity");
      addr[nt] = &All.LogAGNLuminosity;
      id[nt++] = REAL;

      strcpy(tag[nt], "PhotonDesNumNgb");
      addr[nt] = &All.PhotonDesNumNgb;
      id[nt++] = INT;

      strcpy(tag[nt], "PhotonMaxNumNgbDeviation");
      addr[nt] = &All.PhotonMaxNumNgbDeviation;
      id[nt++] = INT;
#endif

#if(defined(MRT_SOURCES) && defined(MRT_STARS)) || (defined(MRT_LOCAL_FEEDBACK))
      strcpy(tag[nt], "SpectrumTablePath");
      addr[nt] = &All.SpectrumTablePath;
      id[nt++] = STRING;
#endif

#if defined(MRT_SOURCES) && defined(MRT_STARS)
      strcpy(tag[nt], "EscapeFraction");
      addr[nt] = &All.EscapeFraction;
      id[nt++] = REAL;
#endif

#endif

#if defined(COOLING) && !defined(SIMPLE_COOLING) && !defined(GRACKLE) && !defined(CHIMES)
      strcpy(tag[nt], "TreecoolFile");
      addr[nt] = &All.TreecoolFile;
      id[nt++] = STRING;
#ifdef UVB_START
      strcpy(tag[nt], "UVBStartRedshift");
      addr[nt] = &All.UVBStartRedshift;
      id[nt++] = REAL;
#endif
#ifdef UVB_SELF_SHIELDING
      strcpy(tag[nt], "SelfShieldingFile");
      addr[nt] = &All.SelfShieldingFile;
      id[nt++] = STRING;
#endif
#ifdef GFM_AGN_RADIATION
      strcpy(tag[nt], "TreecoolFileAGN");
      addr[nt] = &All.TreecoolFileAGN;
      id[nt++] = STRING;
#endif
#ifdef RADCOOL
      strcpy(tag[nt], "TreecoolFileRAD");
      addr[nt] = &All.TreecoolFileRAD;
      id[nt++] = STRING;
#endif
#endif

#ifdef GRACKLE
      strcpy(tag[nt], "GrackleOn");
      addr[nt] = &All.GrackleOn;
      id[nt++] = INT;

      strcpy(tag[nt], "GrackleRadiativeCooling");
      addr[nt] = &All.GrackleRadiativeCooling;
      id[nt++] = INT;

      strcpy(tag[nt], "GrackleMetalCooling");
      addr[nt] = &All.GrackleMetalCooling;
      id[nt++] = INT;

      strcpy(tag[nt], "GrackleUVB");
      addr[nt] = &All.GrackleUVB;
      id[nt++] = INT;

      strcpy(tag[nt], "GrackleDataFile");
      addr[nt] = &All.GrackleDataFile;
      id[nt++] = STRING;
#ifndef METALS
      strcpy(tag[nt], "GrackleInitialMetallicity");
      addr[nt] = &All.GrackleInitialMetallicity;
      id[nt++] = REAL;
#endif
#endif

#ifdef CHIMES
      strcpy(tag[nt], "Chimes_data_path");
      addr[nt] = All.ChimesDataPath;
      id[nt++] = STRING;

      strcpy(tag[nt], "PhotoIon_table_path");
      addr[nt] = ChimesGlobalVars.PhotoIonTablePath;
      id[nt++] = STRING;

      strcpy(tag[nt], "EqAbundance_table_path");
      addr[nt] = ChimesGlobalVars.EqAbundanceTablePath;
      id[nt++] = STRING;

      strcpy(tag[nt], "Thermal_Evolution_On");
      addr[nt] = &All.ChimesThermEvolOn;
      id[nt++] = INT;

      strcpy(tag[nt], "Chemistry_eqm");
      addr[nt] = &All.ChimesForceEqOn;
      id[nt++] = INT;

      strcpy(tag[nt], "Reduction_On");
      addr[nt] = &ChimesGlobalVars.reductionOn;
      id[nt++] = INT;

      strcpy(tag[nt], "StaticMolCooling");
      addr[nt] = &ChimesGlobalVars.StaticMolCooling;
      id[nt++] = INT;

      strcpy(tag[nt], "CellSelfShieldingOn");
      addr[nt] = &ChimesGlobalVars.cellSelfShieldingOn;
      id[nt++] = INT;

#if defined(CHIMES_JEANS_SHIELDING) || defined(CHIMES_SOBOLEV_SHIELDING)
      strcpy(tag[nt], "ShieldingLengthFactor");
      addr[nt] = &All.ChimesShieldingLengthFactor;
      id[nt++] = REAL;

      strcpy(tag[nt], "MaxShieldingLength_kpc");
      addr[nt] = &All.ChimesMaxShieldingLength_kpc;
      id[nt++] = REAL;
#endif

      strcpy(tag[nt], "Reduction_N_Ions_Low");
      addr[nt] = &ChimesGlobalVars.n_ions_low;
      id[nt++] = INT;

      strcpy(tag[nt], "Reduction_N_Ions_Med");
      addr[nt] = &ChimesGlobalVars.n_ions_med;
      id[nt++] = INT;

      strcpy(tag[nt], "Reduction_N_Ions_High");
      addr[nt] = &ChimesGlobalVars.n_ions_high;
      id[nt++] = INT;

      strcpy(tag[nt], "Grain_Temperature");
      addr[nt] = &ChimesGlobalVars.grain_temperature;
      id[nt++] = REAL;

      strcpy(tag[nt], "CrRate");
      addr[nt] = &All.ChimesCrRate;
      id[nt++] = REAL;

      strcpy(tag[nt], "macrostep_tolerance");
      addr[nt] = &ChimesGlobalVars.time_tolerance;
      id[nt++] = REAL;

      strcpy(tag[nt], "min_macrostep");
      addr[nt] = &ChimesGlobalVars.min_subcyclestep;
      id[nt++] = REAL;

      strcpy(tag[nt], "max_mol_temperature");
      addr[nt] = &ChimesGlobalVars.T_mol;
      id[nt++] = REAL;

      strcpy(tag[nt], "IsotropicPhotonDensity");
      addr[nt] = &All.ChimesIsotropicPhotonDensity;
      id[nt++] = REAL;

      strcpy(tag[nt], "TemperatureEqm_Thresh");
      addr[nt] = &ChimesGlobalVars.T_EqThresh;
      id[nt++] = REAL;

      strcpy(tag[nt], "relativeTolerance");
      addr[nt] = &ChimesGlobalVars.relativeTolerance;
      id[nt++] = REAL;

      strcpy(tag[nt], "absoluteTolerance");
      addr[nt] = &ChimesGlobalVars.absoluteTolerance;
      id[nt++] = REAL;

      strcpy(tag[nt], "thermalAbsoluteTolerance");
      addr[nt] = &ChimesGlobalVars.thermalAbsoluteTolerance;
      id[nt++] = REAL;

      strcpy(tag[nt], "scale_metal_tolerances");
      addr[nt] = &ChimesGlobalVars.scale_metal_tolerances;
      id[nt++] = INT;

      strcpy(tag[nt], "IncludeCarbon");
      addr[nt] = &ChimesGlobalVars.element_included[0];
      id[nt++] = INT;

      strcpy(tag[nt], "IncludeNitrogen");
      addr[nt] = &ChimesGlobalVars.element_included[1];
      id[nt++] = INT;

      strcpy(tag[nt], "IncludeOxygen");
      addr[nt] = &ChimesGlobalVars.element_included[2];
      id[nt++] = INT;

      strcpy(tag[nt], "IncludeNeon");
      addr[nt] = &ChimesGlobalVars.element_included[3];
      id[nt++] = INT;

      strcpy(tag[nt], "IncludeMagnesium");
      addr[nt] = &ChimesGlobalVars.element_included[4];
      id[nt++] = INT;

      strcpy(tag[nt], "IncludeSilicon");
      addr[nt] = &ChimesGlobalVars.element_included[5];
      id[nt++] = INT;

      strcpy(tag[nt], "IncludeSulphur");
      addr[nt] = &ChimesGlobalVars.element_included[6];
      id[nt++] = INT;

      strcpy(tag[nt], "IncludeCalcium");
      addr[nt] = &ChimesGlobalVars.element_included[7];
      id[nt++] = INT;

      strcpy(tag[nt], "IncludeIron");
      addr[nt] = &ChimesGlobalVars.element_included[8];
      id[nt++] = INT;
#endif

#if defined(RADCOOL) && defined(COOLING)
      strcpy(tag[nt], "NewStarsOn");
      addr[nt] = &All.NewStarsOn;
      id[nt++] = INT;

      strcpy(tag[nt], "OldStarsOn");
      addr[nt] = &All.OldStarsOn;
      id[nt++] = INT;

#ifdef RADCOOL_HOTHALO
      strcpy(tag[nt], "HotHaloOn");
      addr[nt] = &All.HotHaloOn;
      id[nt++] = INT;
#endif

      strcpy(tag[nt], "SelfShieldingOn");
      addr[nt] = &All.SelfShieldingOn;
      id[nt++] = INT;

      strcpy(tag[nt], "SelfShieldingDensity");
      addr[nt] = &All.SelfShieldingDensity;
      id[nt++] = REAL;
#endif

#ifdef GFM_AGN_RADIATION
      strcpy(tag[nt], "SelfShieldingDensity");
      addr[nt] = &All.SelfShieldingDensity;
      id[nt++] = REAL;

      strcpy(tag[nt], "ObscurationFactor");
      addr[nt] = &All.ObscurationFactor;
      id[nt++] = REAL;

      strcpy(tag[nt], "ObscurationSlope");
      addr[nt] = &All.ObscurationSlope;
      id[nt++] = REAL;
#endif

#if defined(BH_ADIOS_WIND) || defined(BH_BUBBLES)
      strcpy(tag[nt], "RadioFeedbackFactor");
      addr[nt] = &All.RadioFeedbackFactor;
      id[nt++] = REAL;
#endif

#if defined(BH_ADIOS_WIND)
#ifdef BH_ADIOS_ONLY_ABOVE_MINIMUM_DENSITY
      strcpy(tag[nt], "RadioFeedbackMinDensityFactor");
      addr[nt] = &All.RadioFeedbackMinDensityFactor;
      id[nt++] = REAL;
#endif
#ifdef BH_ADIOS_DENS_DEP_EFFICIANCY
      strcpy(tag[nt], "RadioFeedbackFactorPivotMass");
      addr[nt] = &All.RadioFeedbackFactorPivotMass;
      id[nt++] = REAL;

      strcpy(tag[nt], "RadioFeedbackFactorSlope");
      addr[nt] = &All.RadioFeedbackFactorSlope;
      id[nt++] = REAL;

      strcpy(tag[nt], "RadioFeedbackFactorRefDensityFactor");
      addr[nt] = &All.RadioFeedbackFactorRefDensityFactor;
      id[nt++] = REAL;

      strcpy(tag[nt], "RadioFeedbackFactorMaxEfficiency");
      addr[nt] = &All.RadioFeedbackFactorMaxEfficiency;
      id[nt++] = REAL;
#endif

#ifdef BH_ADIOS_RANDOMIZED
      strcpy(tag[nt], "RadioFeedbackReiorientationFactor");
      addr[nt] = &All.RadioFeedbackReiorientationFactor;
      id[nt++] = REAL;
#endif

#if defined(BH_ADIOS_WIND_WITH_QUASARTHRESHOLD)
      strcpy(tag[nt], "QuasarThreshold");
      addr[nt] = &All.QuasarThreshold;
      id[nt++] = REAL;
#endif
#endif

#ifdef BH_BUBBLES
      strcpy(tag[nt], "BubbleDistance");
      addr[nt] = &All.BubbleDistance;
      id[nt++] = REAL;

      strcpy(tag[nt], "BubbleRadius");
      addr[nt] = &All.BubbleRadius;
      id[nt++] = REAL;

      strcpy(tag[nt], "BubbleEnergy");
      addr[nt] = &All.BubbleEnergy;
      id[nt++] = REAL;

      strcpy(tag[nt], "BlackHoleRadioTriggeringFactor");
      addr[nt] = &All.BlackHoleRadioTriggeringFactor;
      id[nt++] = REAL;

      strcpy(tag[nt], "DefaultICMDensity");
      addr[nt] = &All.DefaultICMDensity;
      id[nt++] = REAL;

#ifdef UNIFIED_FEEDBACK
      strcpy(tag[nt], "RadioThreshold");
      addr[nt] = &All.RadioThreshold;
      id[nt++] = REAL;
#endif
#ifdef BH_MAGNETIC_BUBBLES
      strcpy(tag[nt], "MagneticEnergyFraction");
      addr[nt] = &All.MagneticEnergyFraction;
      id[nt++] = REAL;
#endif
#endif

#ifdef GFM_UVB_CORRECTIONS
      strcpy(tag[nt], "UV_HeII_threshold");
      addr[nt] = &All.UV_HeII_threshold;
      id[nt++] = REAL;

      strcpy(tag[nt], "UV_HeII_alpha");
      addr[nt] = &All.UV_HeII_alpha;
      id[nt++] = REAL;

      strcpy(tag[nt], "UV_HeII_beta");
      addr[nt] = &All.UV_HeII_alpha;
      id[nt++] = REAL;

      strcpy(tag[nt], "UV_HeatBoost");
      addr[nt] = &All.UV_HeatBoost;
      id[nt++] = REAL;
#endif

#ifdef GFM_WIND_ENERGY_METAL_DEPENDENCE
      strcpy(tag[nt], "WindEnergyReductionFactor");
      addr[nt] = &All.WindEnergyReductionFactor;
      id[nt++] = REAL;

      strcpy(tag[nt], "WindEnergyReductionMetallicity");
      addr[nt] = &All.WindEnergyReductionMetallicity;
      id[nt++] = REAL;

      strcpy(tag[nt], "WindEnergyReductionExponent");
      addr[nt] = &All.WindEnergyReductionExponent;
      id[nt++] = REAL;
#endif

#ifdef GFM_STELLAR_EVOLUTION
      strcpy(tag[nt], "DesNumNgbEnrichment");
      addr[nt] = &All.DesNumNgbEnrichment;
      id[nt++] = INT;

      strcpy(tag[nt], "MaxNumNgbDeviationEnrichment");
      addr[nt] = &All.MaxNumNgbDeviationEnrichment;
      id[nt++] = REAL;

      strcpy(tag[nt], "SNII_MinMass_Msun");
      addr[nt] = &All.SNII_MinMass_Msun;
      id[nt++] = REAL;

      strcpy(tag[nt], "SNII_MaxMass_Msun");
      addr[nt] = &All.SNII_MaxMass_Msun;
      id[nt++] = REAL;

      strcpy(tag[nt], "IMF_MinMass_Msun");
      addr[nt] = &All.IMF_MinMass_Msun;
      id[nt++] = REAL;

      strcpy(tag[nt], "IMF_MaxMass_Msun");
      addr[nt] = &All.IMF_MaxMass_Msun;
      id[nt++] = REAL;

      strcpy(tag[nt], "AGB_MassTransferOn");
      addr[nt] = &All.AGB_MassTransferOn;
      id[nt++] = INT;

      strcpy(tag[nt], "SNIa_Rate_Norm");
      addr[nt] = &All.SNIa_Rate_Norm;
      id[nt++] = REAL;

      strcpy(tag[nt], "SNIa_Rate_TAU");
      addr[nt] = &All.SNIa_Rate_TAU;
      id[nt++] = REAL;

      strcpy(tag[nt], "SNIa_MassTransferOn");
      addr[nt] = &All.SNIa_MassTransferOn;
      id[nt++] = INT;

      strcpy(tag[nt], "SNII_MassTransferOn");
      addr[nt] = &All.SNII_MassTransferOn;
      id[nt++] = INT;

#ifdef GFM_RPROCESS
      strcpy(tag[nt], "NSNS_MassTransferOn");
      addr[nt] = &All.NSNS_MassTransferOn;
      id[nt++] = INT;

      strcpy(tag[nt], "NSNS_Rate_TAU");
      addr[nt] = &All.NSNS_Rate_TAU;
      id[nt++] = REAL;

      strcpy(tag[nt], "NSNS_MassPerEvent");
      addr[nt] = &All.NSNS_MassPerEvent;
      id[nt++] = REAL;

      strcpy(tag[nt], "NSNS_per_SNIa");
      addr[nt] = &All.NSNS_per_SNIa;
      id[nt++] = REAL;
#endif

#ifdef GFM_RPROCESS_CHANNELS
      for(i = 0; i < GFM_RPROCESS_NSNS; i++)
        {
          sprintf(buf, "NSNS_MassPerEvent%d", i);
          strcpy(tag[nt], buf);
          addr[nt] = &All.rp[i].NSNS_MassPerEvent;
          id[nt++] = REAL;

          sprintf(buf, "NSNS_RateNorm%d", i);
          strcpy(tag[nt], buf);
          addr[nt] = &All.rp[i].NSNS_RateNorm;
          id[nt++] = REAL;

          sprintf(buf, "NSNS_RateTAU%d", i);
          strcpy(tag[nt], buf);
          addr[nt] = &All.rp[i].NSNS_RateTAU;
          id[nt++] = REAL;

          sprintf(buf, "NSNS_PowerlawIndex%d", i);
          strcpy(tag[nt], buf);
          addr[nt] = &All.rp[i].NSNS_PowerlawIndex;
          id[nt++] = REAL;

#ifdef GFM_RPROCESS_CHANNELS_NS_KICKS
          sprintf(buf, "NSNS_KickAverage%d", i);
          strcpy(tag[nt], buf);
          addr[nt] = &All.rp[i].NSNS_KickAverage;
          id[nt++] = REAL;

          sprintf(buf, "NSNS_KickSigma%d", i);
          strcpy(tag[nt], buf);
          addr[nt] = &All.rp[i].NSNS_KickSigma;
          id[nt++] = REAL;
#endif
        }

#if GFM_RPROCESS_NSNS < GFM_RPROCESS_CHANNELS
      for(i = GFM_RPROCESS_NSNS; i < GFM_RPROCESS_CHANNELS; i++)
        {
          sprintf(buf, "RPSN_MassPerEvent%d", i - GFM_RPROCESS_NSNS);
          strcpy(tag[nt], buf);
          addr[nt] = &All.rpSN[i - GFM_RPROCESS_NSNS].RPSN_MassPerEvent;
          id[nt++] = REAL;

          sprintf(buf, "RPSN_FractionPerSN%d", i - GFM_RPROCESS_NSNS);
          strcpy(tag[nt], buf);
          addr[nt] = &All.rpSN[i - GFM_RPROCESS_NSNS].RPSN_FractionPerSN;
          id[nt++] = REAL;
        }
#endif
#endif

      strcpy(tag[nt], "YieldTablePath");
      addr[nt] = &All.YieldTablePath;
      id[nt++] = STRING;
#endif

#ifdef GFM_DUST
      strcpy(tag[nt], "AGB_Dust_Delta_C");
      addr[nt] = &All.AGB_Dust_Delta_C;
      id[nt++] = REAL;

      strcpy(tag[nt], "AGB_Dust_Delta_Metal");
      addr[nt] = &All.AGB_Dust_Delta_Metal;
      id[nt++] = REAL;

      strcpy(tag[nt], "SN_Dust_Delta_C");
      addr[nt] = &All.SN_Dust_Delta_C;
      id[nt++] = REAL;

      strcpy(tag[nt], "SN_Dust_Delta_Metal");
      addr[nt] = &All.SN_Dust_Delta_Metal;
      id[nt++] = REAL;

      strcpy(tag[nt], "Dust_Growth_Tau");
      addr[nt] = &All.Dust_Growth_Tau;
      id[nt++] = REAL;

#ifndef GFM_DUST_MRN
      strcpy(tag[nt], "DustSingleGrainSize");
      addr[nt] = &All.DustSingleGrainSize;
      id[nt++] = REAL;
#endif

#ifdef GFM_DUST_SPUTTERING
#if(GFM_DUST_SPUTTERING == 1)
      strcpy(tag[nt], "Dust_Sputter_Tau_Fac");
      addr[nt] = &All.Dust_Sputter_Tau_Fac;
      id[nt++] = REAL;
#endif
#endif

#if(GFM_DUST_DESTMODE == 1)
      strcpy(tag[nt], "Dust_Destruction_Tau");
      addr[nt] = &All.Dust_Destruction_Tau;
      id[nt++] = REAL;
#endif
#endif

#ifdef GFM_PREENRICH
#ifndef CHIMES_PREENRICH_AT_START
      strcpy(tag[nt], "PreEnrichTime");
      addr[nt] = &All.PreEnrichTime;
      id[nt++] = INT;
#endif

      strcpy(tag[nt], "PreEnrichAbundanceFile");
      addr[nt] = &All.PreEnrichAbundanceFile;
      id[nt++] = STRING;
#endif

#ifdef GFM_COOLING_METAL
      strcpy(tag[nt], "CoolingTablePath");
      addr[nt] = &All.CoolingTablePath;
      id[nt++] = STRING;

      strcpy(tag[nt], "MinMetalTemp");
      addr[nt] = &All.MinMetalTemp;
      id[nt++] = REAL;
#endif

#ifdef GFM_STELLAR_PHOTOMETRICS
      strcpy(tag[nt], "PhotometricsTablePath");
      addr[nt] = &All.PhotometricsTablePath;
      id[nt++] = STRING;
#endif

#if defined(REFINEMENT) || defined(GFM_STELLAR_EVOLUTION) || defined(BLACK_HOLES)
      strcpy(tag[nt], "ReferenceGasPartMass");
      addr[nt] = &All.ReferenceGasPartMass;
      id[nt++] = REAL;
#endif

#if defined(REFINEMENT) && !defined(DG)
      strcpy(tag[nt], "TargetGasMassFactor");
      addr[nt] = &All.TargetGasMassFactor;
      id[nt++] = REAL;

      strcpy(tag[nt], "RefinementCriterion");
      addr[nt] = &All.RefinementCriterion;
      id[nt++] = INT;

      strcpy(tag[nt], "DerefinementCriterion");
      addr[nt] = &All.DerefinementCriterion;
      id[nt++] = INT;
#endif

#ifdef GMC_REFINEMENT
      strcpy(tag[nt], "GMCRefCellsPerJeansLength");
      addr[nt] = &All.GMCRefCellsPerJeansLength;
      id[nt++] = REAL;

      strcpy(tag[nt], "GMCRefMinDensity");
      addr[nt] = &All.GMCRefMinDensity;
      id[nt++] = REAL;

      strcpy(tag[nt], "GMCRefMaxDensity");
      addr[nt] = &All.GMCRefMaxDensity;
      id[nt++] = REAL;

      strcpy(tag[nt], "GMCDerefMinDensity");
      addr[nt] = &All.GMCDerefMinDensity;
      id[nt++] = REAL;
#endif

#ifdef STICKYFLAGS
#if(REFLECTIVE_X == 2) || (REFLECTIVE_Y == 2) || (REFLECTIVE_Z == 2)
      strcpy(tag[nt], "StickyLayerMaxDist");
      addr[nt] = &All.StickyLayerMaxDist;
      id[nt++] = REAL;
#endif
#endif

#ifdef GFM_STELLAR_FEEDBACK
      strcpy(tag[nt], "EnergyPerSNIa");
      addr[nt] = &All.EnergyPerSNIa;
      id[nt++] = REAL;

      strcpy(tag[nt], "AGBWindVelocity");
      addr[nt] = &All.AGBWindVelocity;
      id[nt++] = REAL;
#endif

#ifdef USE_SFR

#ifdef QUICK_LYALPHA_LATETIMEONLY
      strcpy(tag[nt], "TimeOfCoolingStart");
      addr[nt] = &All.TimeOfCoolingStart;
      id[nt++] = REAL;
#endif

#ifdef STEEPER_SFR_FOR_STARBURST
      strcpy(tag[nt], "StarburstPowerLawIndex");
      addr[nt] = &All.StarburstPowerLawIndex;
      id[nt++] = REAL;
#endif

#if !defined(FM_SFR) && !defined(LOCAL_FEEDBACK) && !defined(SFR_MCS)
      strcpy(tag[nt], "CritOverDensity");
      addr[nt] = &All.CritOverDensity;
      id[nt++] = REAL;

      strcpy(tag[nt], "TemperatureThresh");
      addr[nt] = &All.TemperatureThresh;
      id[nt++] = REAL;

      strcpy(tag[nt], "CritPhysDensity");
      addr[nt] = &All.CritPhysDensity;
      id[nt++] = REAL;

#if !defined(GFM_STELLAR_EVOLUTION) || \
    (defined(GFM_STELLAR_EVOLUTION) && defined(GFM_VARIABLE_IMF))  // the latter case is needed for the eEOS (and only there)
      strcpy(tag[nt], "FactorSN");
      addr[nt] = &All.FactorSN;
      id[nt++] = REAL;
#endif

      strcpy(tag[nt], "FactorEVP");
      addr[nt] = &All.FactorEVP;
      id[nt++] = REAL;

#ifdef SOFTEREQS
      strcpy(tag[nt], "FactorForSofterEQS");
      addr[nt] = &All.FactorForSofterEQS;
      id[nt++] = REAL;

      strcpy(tag[nt], "TempForSofterEQS");
      addr[nt] = &All.TempForSofterEQS;
      id[nt++] = REAL;
#endif

#ifdef MODIFIED_EOS
      strcpy(tag[nt], "FactorDensThresh");
      addr[nt] = &All.FactorDensThresh;
      id[nt++] = REAL;

      strcpy(tag[nt], "FactorUthermAtThresh");
      addr[nt] = &All.FactorUthermAtThresh;
      id[nt++] = REAL;

      strcpy(tag[nt], "FactorUthermJoin");
      addr[nt] = &All.FactorUthermJoin;
      id[nt++] = REAL;
#endif

      strcpy(tag[nt], "TempSupernova");
      addr[nt] = &All.TempSupernova;
      id[nt++] = REAL;

      strcpy(tag[nt], "TempClouds");
      addr[nt] = &All.TempClouds;
      id[nt++] = REAL;

      strcpy(tag[nt], "MaxSfrTimescale");
      addr[nt] = &All.MaxSfrTimescale;
      id[nt++] = REAL;

// #endif
#endif

#ifdef GFM_WINDS_LOCAL
      strcpy(tag[nt], "VariableWindVelFactor");
      addr[nt] = &All.VariableWindVelFactor;
      id[nt++] = REAL;

      strcpy(tag[nt], "WindEnergyIn1e51erg");
      addr[nt] = &All.WindEnergyIn1e51erg;
      id[nt++] = REAL;

      strcpy(tag[nt], "WindFreeTravelMaxTimeFactor");
      addr[nt] = &All.WindFreeTravelMaxTimeFactor;
      id[nt++] = REAL;

      strcpy(tag[nt], "WindFreeTravelDensFac");
      addr[nt] = &All.WindFreeTravelDensFac;
      id[nt++] = REAL;

      strcpy(tag[nt], "MinWindVel");
      addr[nt] = &All.MinWindVel;
      id[nt++] = REAL;
#endif

#if defined(GFM_CONST_IMF) && (GFM_CONST_IMF == 1)
      strcpy(tag[nt], "IMFslope");
      addr[nt] = &All.IMFslope;
      id[nt++] = REAL;
#endif

#ifdef GFM_WINDS_STRIPPING
      strcpy(tag[nt], "WindDumpFactor");
      addr[nt] = &All.WindDumpFactor;
      id[nt++] = REAL;
#endif

#ifdef GFM_WINDS
#ifndef GFM_WINDS_VARIABLE
      strcpy(tag[nt], "WindEfficiency");
      addr[nt] = &All.WindEfficiency;
      id[nt++] = REAL;
#endif

#ifdef GFM_STELLAR_EVOLUTION
      strcpy(tag[nt], "WindEnergyIn1e51erg");
      addr[nt] = &All.WindEnergyIn1e51erg;
      id[nt++] = REAL;
#else
      strcpy(tag[nt], "WindEnergyFraction");
      addr[nt] = &All.WindEnergyFraction;
      id[nt++] = REAL;
#endif

      strcpy(tag[nt], "WindFreeTravelMaxTimeFactor");
      addr[nt] = &All.WindFreeTravelMaxTimeFactor;
      id[nt++] = REAL;

      strcpy(tag[nt], "WindFreeTravelDensFac");
      addr[nt] = &All.WindFreeTravelDensFac;
      id[nt++] = REAL;

#ifdef GFM_WINDS_THERMAL
      strcpy(tag[nt], "ThermalWindFactor");
      addr[nt] = &All.ThermalWindFactor;
      id[nt++] = REAL;
#endif

#ifdef GFM_WINDS_THERMAL_NEWDEF
      strcpy(tag[nt], "ThermalWindFraction");
      addr[nt] = &All.ThermalWindFraction;
      id[nt++] = REAL;
#endif

#ifdef GFM_WINDS_VARIABLE
      strcpy(tag[nt], "VariableWindVelFactor");
      addr[nt] = &All.VariableWindVelFactor;
      id[nt++] = REAL;

      strcpy(tag[nt], "VariableWindSpecMomentum");
      addr[nt] = &All.VariableWindSpecMomentum;
      id[nt++] = REAL;

      strcpy(tag[nt], "MinWindVel");
      addr[nt] = &All.MinWindVel;
      id[nt++] = REAL;

#ifdef GFM_WINDS_MASSSCALING
      strcpy(tag[nt], "VariableWindMassScale");
      addr[nt] = &All.VariableWindMassScale;
      id[nt++] = REAL;
#endif

#ifdef GFM_WINDS_HUBBLESCALING
      strcpy(tag[nt], "WindSuppressionRedshift");
      addr[nt] = &All.WindSuppressionRedshift;
      id[nt++] = REAL;
#endif

#if(GFM_WINDS_VARIABLE == 0)
      strcpy(tag[nt], "HaloConcentrationNorm");
      addr[nt] = &All.HaloConcentrationNorm;
      id[nt++] = REAL;

      strcpy(tag[nt], "HaloConcentrationSlope");
      addr[nt] = &All.HaloConcentrationSlope;
      id[nt++] = REAL;
#endif
#endif
#endif
#endif

#ifdef FM_SFR
      strcpy(tag[nt], "CritOverDensity");
      addr[nt] = &All.CritOverDensity;
      id[nt++] = REAL;

      strcpy(tag[nt], "DensThreshold");
      addr[nt] = &All.DensThreshold;
      id[nt++] = REAL;

      strcpy(tag[nt], "SfrEfficiency");
      addr[nt] = &All.SfrEfficiency;
      id[nt++] = REAL;

#ifdef USE_POLYTROPIC_EQSTATE
      strcpy(tag[nt], "TemperatureAtThreshold");
      addr[nt] = &All.UthermThreshold;
      id[nt++] = REAL;
#endif
#endif

#ifdef FM_STAR_FEEDBACK
      strcpy(tag[nt], "FeedbackEfficiency");
      addr[nt] = &All.FeedbackEfficiency;
      id[nt++] = REAL;

      strcpy(tag[nt], "EtaKineticEnergy");
      addr[nt] = &All.EtaKineticEnergy;
      id[nt++] = REAL;

      strcpy(tag[nt], "PhotoionizationGasTemp");
      addr[nt] = &All.PhotoionizationGasTemp;
      id[nt++] = REAL;

#if !defined(DELAYED_COOLING) && !defined(DELAYED_COOLING_TURB) && !defined(NON_STOCHASTIC_MOMENTUM_FEEDBACK) && \
    !defined(DIRECT_MOMENTUM_INJECTION_FEEDBACK)
      strcpy(tag[nt], "FeedbackVelocity");
      addr[nt] = &All.FeedbackVelocity;
      id[nt++] = REAL;
#endif

#ifdef DELAYED_COOLING_TURB
      strcpy(tag[nt], "DissipationTime");
      addr[nt] = &All.DissipationTime;
      id[nt++] = REAL;

      strcpy(tag[nt], "SigmaThreshold");
      addr[nt] = &All.SigmaThreshold;
      id[nt++] = REAL;
#endif
#ifdef INSTANTANEOUS_DEPOSITION
      strcpy(tag[nt], "MinFeedAge");
      addr[nt] = &All.MinFeedAge;
      id[nt++] = REAL;

      strcpy(tag[nt], "FeedbackInjectionEvents");
      addr[nt] = &All.FeedbackInjectionEvents;
      id[nt++] = INT;
#endif
#ifdef DIRECT_MOMENTUM_INJECTION_FEEDBACK
      strcpy(tag[nt], "MomentumFactor");
      addr[nt] = &All.MomentumFactor;
      id[nt++] = REAL;
#endif
#endif

#ifdef FM_RADIATION_FEEDBACK
      strcpy(tag[nt], "DustOppacityRadiationFeedback");
      addr[nt] = &All.DustOppacityRadiationFeedback;
      id[nt++] = REAL;

      strcpy(tag[nt], "InputTimeHeatRadiationFeedback");
      addr[nt] = &All.InputTimeHeatRadiationFeedback;
      id[nt++] = REAL;

      strcpy(tag[nt], "InputTimeMomRadiationFeedback");
      addr[nt] = &All.InputTimeMomRadiationFeedback;
      id[nt++] = REAL;

      strcpy(tag[nt], "LumToMassRatioRadiationFeedback");
      addr[nt] = &All.LumToMassRatioRadiationFeedback;
      id[nt++] = REAL;

      strcpy(tag[nt], "RadiationFeedbackAvgPhotonEnergyineV");
      addr[nt] = &All.RadiationFeedbackAvgPhotonEnergyineV;
      id[nt++] = REAL;
#endif
#ifdef FM_EARLY_STAR_FEEDBACK
      strcpy(tag[nt], "EarlyEtaKineticEnergy");
      addr[nt] = &All.EarlyEtaKineticEnergy;
      id[nt++] = REAL;

      strcpy(tag[nt], "EarlyFeedbackEfficiency");
      addr[nt] = &All.EarlyFeedbackEfficiency;
      id[nt++] = REAL;

/*
      strcpy(tag[nt], "EarlyFeedbackVelocity");
      addr[nt] = &All.EarlyFeedbackVelocity;
      id[nt++] = REAL;
*/
#endif

#ifdef SFR_MCS
      strcpy(tag[nt], "DensThreshold");
      addr[nt] = &All.DensThreshold;
      id[nt++] = REAL;

      strcpy(tag[nt], "SfrEfficiency");
      addr[nt] = &All.SfrEfficiency;
      id[nt++] = REAL;

      strcpy(tag[nt], "CritOverDensity");
      addr[nt] = &All.CritOverDensity;
      id[nt++] = REAL;

#ifdef SN_MCS
      strcpy(tag[nt], "MaxSNEvalTimestep");
      addr[nt] = &All.MaxSNEvalTimestep;
      id[nt++] = REAL;

      strcpy(tag[nt], "MaxSNStarMassFac");
      addr[nt] = &All.MaxSNStarMassFac;
      id[nt++] = REAL;

      strcpy(tag[nt], "HostCellFeedbackFraction");
      addr[nt] = &All.HostCellFeedbackFraction;
      id[nt++] = REAL;

      strcpy(tag[nt], "SB99TablesPath");
      addr[nt] = &All.SB99TablesPath;
      id[nt++] = STRING;

#ifdef SB99_FIXED_Z
      strcpy(tag[nt], "SB99_metallicity_name");
      addr[nt] = &All.SB99_metallicity_name;
      id[nt++] = STRING;
#endif

      strcpy(tag[nt], "SNKineticRatio");
      addr[nt] = &All.SNKineticRatio;
      id[nt++] = REAL;

      strcpy(tag[nt], "SNMassReturn");
      addr[nt] = &All.SNMassReturn;
      id[nt++] = REAL;
#ifdef METALS
      strcpy(tag[nt], "SNEjectaMetallicity");
      addr[nt] = &All.SNEjectaMetallicity;
      id[nt++] = REAL;
#endif
#endif
#endif

#ifdef DARKENERGY
#ifndef TIMEDEPDE
      strcpy(tag[nt], "DarkEnergyParam");
      addr[nt] = &All.DarkEnergyParam;
      id[nt++] = REAL;
#endif
#endif

#ifdef DARKENERGY
#ifdef TIMEDEPDE
      strcpy(tag[nt], "DarkEnergyFile");
      addr[nt] = All.DarkEnergyFile;
      id[nt++] = STRING;
#endif
#endif

#ifdef RELAXOBJECT
      strcpy(tag[nt], "RelaxBaseFac");
      addr[nt] = &All.RelaxBaseFac;
      id[nt++] = REAL;
#endif
#ifdef RELAXOBJECT_COOLING
      strcpy(tag[nt], "RelaxTemperature");
      addr[nt] = &All.RelaxTemperature;
      id[nt++] = REAL;
#endif
#ifdef RELAXOBJECT_COOLING2
      strcpy(tag[nt], "TempCore");
      addr[nt] = &All.TempCore;
      id[nt++] = REAL;

      strcpy(tag[nt], "TempShell");
      addr[nt] = &All.TempShell;
      id[nt++] = REAL;

      strcpy(tag[nt], "MassCore");
      addr[nt] = &All.MassCore;
      id[nt++] = REAL;

      strcpy(tag[nt], "ShellBaseRadius");
      addr[nt] = &All.ShellBaseRadius;
      id[nt++] = REAL;

      strcpy(tag[nt], "TempMin");
      addr[nt] = &All.TempMin;
      id[nt++] = REAL;

      strcpy(tag[nt], "RTempMin");
      addr[nt] = &All.RTempMin;
      id[nt++] = REAL;

      strcpy(tag[nt], "BackgroundTemp");
      addr[nt] = &All.BackgroundTemp;
      id[nt++] = REAL;
#endif

#ifdef INSPIRAL
      strcpy(tag[nt], "InspiralVelocity");
      addr[nt] = &All.InspiralVelocity;
      id[nt++] = REAL;
#endif

#ifdef EOS_DEGENERATE
      strcpy(tag[nt], "EosTable");
      addr[nt] = All.EosTable;
      id[nt++] = STRING;

      strcpy(tag[nt], "EosSpecies");
      addr[nt] = All.EosSpecies;
      id[nt++] = STRING;
#endif

#ifdef NUCLEAR_NETWORK
      strcpy(tag[nt], "NetworkRates");
      addr[nt] = All.NetworkRates;
      id[nt++] = STRING;

      strcpy(tag[nt], "NetworkPartFunc");
      addr[nt] = All.NetworkPartFunc;
      id[nt++] = STRING;

      strcpy(tag[nt], "NetworkMasses");
      addr[nt] = All.NetworkMasses;
      id[nt++] = STRING;

      strcpy(tag[nt], "NetworkWeakrates");
      addr[nt] = All.NetworkWeakrates;
      id[nt++] = STRING;

      strcpy(tag[nt], "NetworkTempThreshold");
      addr[nt] = &All.NetworkTempThreshold;
      id[nt++] = REAL;

#ifdef NUCLEAR_NETWORK_LIMIT_COMPOSITION_CHANGE
      strcpy(tag[nt], "NuclearNetworkMaxCompositionChange");
      addr[nt] = &All.NuclearNetworkMaxCompositionChange;
      id[nt++] = REAL;
#else
      strcpy(tag[nt], "NuclearNetworkMaxTempChange");
      addr[nt] = &All.NuclearNetworkMaxTempChange;
      id[nt++] = REAL;
#endif

#ifdef NUCLEAR_NETWORK_TIMESTEP_LIMITER
      strcpy(tag[nt], "NuclearNetworkMaxEnergyDiff");
      addr[nt] = &All.NuclearNetworkMaxEnergyDiff;
      id[nt++] = REAL;
#endif
#endif

#ifdef NETWORK_NSE
      strcpy(tag[nt], "NetworkNSEThreshold");
      addr[nt] = &All.NetworkNSEThreshold;
      id[nt++] = REAL;
#endif

#ifdef EOS_OPAL
      strcpy(tag[nt], "EosOpalTable");
      addr[nt] = All.EosOpalTable;
      id[nt++] = STRING;
#endif

#ifdef TRACER_TRAJECTORY
#ifdef TRACER_TRAJECTORY_GENERATE
      strcpy(tag[nt], "NumberOfTracersToGenerate");
      addr[nt] = &All.NumberOfTracersToGenerate;
      id[nt++] = INT;
#else
      strcpy(tag[nt], "TracerInitFile");
      addr[nt] = All.TracerInitFile;
      id[nt++] = STRING;
#endif

      strcpy(tag[nt], "TracerOutputFile");
      addr[nt] = All.TracerOutputFile;
      id[nt++] = STRING;

      strcpy(tag[nt], "TracerOutputConfFile");
      addr[nt] = All.TracerOutputConfFile;
      id[nt++] = STRING;
#endif

#ifdef MHD_CT
      strcpy(tag[nt], "CTMeanBx");
      addr[nt] = &All.CT_mean_Bx;
      id[nt++] = REAL;

      strcpy(tag[nt], "CTMeanBy");
      addr[nt] = &All.CT_mean_By;
      id[nt++] = REAL;

      strcpy(tag[nt], "CTMeanBz");
      addr[nt] = &All.CT_mean_Bz;
      id[nt++] = REAL;
#endif

#ifdef MHD_SEEDFIELD
      strcpy(tag[nt], "MHDSeedDir");
      addr[nt] = &All.B_dir;
      id[nt++] = INT;

      strcpy(tag[nt], "MHDSeedValue");
      addr[nt] = &All.B_value;
      id[nt++] = REAL;
#endif

#ifdef MHD_SEEDPSPEC
      strcpy(tag[nt], "MHDSeedPSpecSlope");
      addr[nt] = &All.B_pspec_slope;
      id[nt++] = REAL;

      strcpy(tag[nt], "MHDSeedPSpecAmpl");
      addr[nt] = &All.B_pspec_ampl;
      id[nt++] = REAL;

      strcpy(tag[nt], "MHDSeedPSpecKCut");
      addr[nt] = &All.B_pspec_kcut;
      id[nt++] = REAL;

      strcpy(tag[nt], "MHDSeedPSpecHelical");
      addr[nt] = &All.B_pspec_helical;
      id[nt++] = REAL;
#endif

#ifdef DMPIC
      strcpy(tag[nt], "PicList");
      addr[nt] = All.PicList;
      id[nt++] = STRING;

      strcpy(tag[nt], "PicDataDir");
      addr[nt] = All.PicDataDir;
      id[nt++] = STRING;

      strcpy(tag[nt], "PicDataName");
      addr[nt] = All.PicDataName;
      id[nt++] = STRING;

      strcpy(tag[nt], "PixelsX");
      addr[nt] = &All.PixelsX;
      id[nt++] = INT;

      strcpy(tag[nt], "PixelsY");
      addr[nt] = &All.PixelsY;
      id[nt++] = INT;
#endif

#ifdef TGSET
      strcpy(tag[nt], "MaxTimeBinDiff");
      addr[nt] = &TGD.MaxTimeBinDiff;
      id[nt++] = INT;

      strcpy(tag[nt], "SnapFreeFallTimeFac");
      addr[nt] = &TGD.SnapFreeFallTimeFac;
      id[nt++] = REAL;

      strcpy(tag[nt], "NHInit");
      addr[nt] = &TGD.NHInit;
      id[nt++] = REAL;

      strcpy(tag[nt], "NHTerm");
      addr[nt] = &TGD.NHTerm;
      id[nt++] = REAL;

      strcpy(tag[nt], "StreamVel");
      addr[nt] = &TGD.StreamVel;
      id[nt++] = REAL;

      strcpy(tag[nt], "JeansTemp");
      addr[nt] = &TGD.JeansTemp;
      id[nt++] = REAL;

      strcpy(tag[nt], "JeansDensLow");
      addr[nt] = &TGD.JeansDensLow;
      id[nt++] = REAL;

      strcpy(tag[nt], "JeansDensHigh");
      addr[nt] = &TGD.JeansDensHigh;
      id[nt++] = REAL;

      strcpy(tag[nt], "JeansNumberLow");
      addr[nt] = &TGD.JeansNumberLow;
      id[nt++] = REAL;

      strcpy(tag[nt], "JeansNumberHigh");
      addr[nt] = &TGD.JeansNumberHigh;
      id[nt++] = REAL;
#endif

#ifdef TGCHEM
      strcpy(tag[nt], "ChemMode");
      addr[nt] = &TGCD.ChemMode;
      id[nt++] = INT;

      strcpy(tag[nt], "ChemIOMode");
      addr[nt] = &TGCD.ChemIOMode;
      id[nt++] = INT;

      strcpy(tag[nt], "ChemH2Mode");
      addr[nt] = &TGCD.ChemH2Mode;
      id[nt++] = INT;

      strcpy(tag[nt], "ChemInitAbH2");
      addr[nt] = &TGCD.ChemInitAbH2;
      id[nt++] = REAL;

      strcpy(tag[nt], "ChemInitAbHII");
      addr[nt] = &TGCD.ChemInitAbHII;
      id[nt++] = REAL;

      strcpy(tag[nt], "ChemJ21");
      addr[nt] = &TGCD.ChemJ21;
      id[nt++] = REAL;
#endif

#ifdef HEALRAY
      strcpy(tag[nt], "RaySourceFile");
      addr[nt] = &HRD.RaySourceFile;
      id[nt++] = STRING;

      strcpy(tag[nt], "RayMode");
      addr[nt] = &HRD.RayMode;
      id[nt++] = INT;

      strcpy(tag[nt], "RayIOMode");
      addr[nt] = &HRD.RayIOMode;
      id[nt++] = INT;

      strcpy(tag[nt], "RayNumThreads");
      addr[nt] = &HRD.RayNumThreads;
      id[nt++] = INT;

      strcpy(tag[nt], "RayDebugMode");
      addr[nt] = &HRD.RayDebugMode;
      id[nt++] = INT;

      strcpy(tag[nt], "RayTimeFac");
      addr[nt] = &HRD.RayTimeFac;
      id[nt++] = INT;

      strcpy(tag[nt], "RayMultiFreqMode");
      addr[nt] = &HRD.RayMultiFreqMode;
      id[nt++] = INT;

      strcpy(tag[nt], "RayNumSources");
      addr[nt] = &HRD.RayNumSources;
      id[nt++] = INT;

      strcpy(tag[nt], "RayNumTrace");
      addr[nt] = &HRD.RayNumTrace;
      id[nt++] = INT;

      strcpy(tag[nt], "RayLvlInit");
      addr[nt] = &HRD.RayLvlInit;
      id[nt++] = INT;

      strcpy(tag[nt], "RayNumLines");
      addr[nt] = &HRD.RayNumLines;
      id[nt++] = INT;

      strcpy(tag[nt], "RayNumFreq");
      addr[nt] = &HRD.RayNumFreq;
      id[nt++] = INT;

      strcpy(tag[nt], "RaySplitFac");
      addr[nt] = &HRD.RaySplitFac;
      id[nt++] = REAL;

      strcpy(tag[nt], "RayPartFrac");
      addr[nt] = &HRD.RayPartFrac;
      id[nt++] = REAL;

      strcpy(tag[nt], "RayMinEgyFrac");
      addr[nt] = &HRD.RayMinEgyFrac;
      id[nt++] = REAL;
#endif

#ifdef SGCHEM
      strcpy(tag[nt], "SGChemConstInitAbundances");
      addr[nt] = &All.SGChemConstInitAbundances;
      id[nt++] = INT;

      strcpy(tag[nt], "SGChemInitH2Abund");
      addr[nt] = &All.SGChemInitH2Abund;
      id[nt++] = REAL;

      mpi_printf("READING chemical variable parameters\n");

      strcpy(tag[nt], "SGChemInitHPAbund");
      addr[nt] = &All.SGChemInitHPAbund;
      id[nt++] = REAL;

#if CHEMISTRYNETWORK == 1
      strcpy(tag[nt], "SGChemInitDIIAbund");
      addr[nt] = &All.SGChemInitDIIAbund;
      id[nt++] = REAL;

      strcpy(tag[nt], "SGChemInitHDAbund");
      addr[nt] = &All.SGChemInitHDAbund;
      id[nt++] = REAL;

      strcpy(tag[nt], "SGChemInitHeIIIAbund");
      addr[nt] = &All.SGChemInitHeIIIAbund;
      id[nt++] = REAL;
#endif

      strcpy(tag[nt], "SGChemInitCPAbund");
      addr[nt] = &All.SGChemInitCPAbund;
      id[nt++] = REAL;

      strcpy(tag[nt], "SGChemInitCOAbund");
      addr[nt] = &All.SGChemInitCOAbund;
      id[nt++] = REAL;

      strcpy(tag[nt], "SGChemInitCHxAbund");
      addr[nt] = &All.SGChemInitCHxAbund;
      id[nt++] = REAL;

      strcpy(tag[nt], "SGChemInitOHxAbund");
      addr[nt] = &All.SGChemInitOHxAbund;
      id[nt++] = REAL;

      strcpy(tag[nt], "SGChemInitHCOPAbund");
      addr[nt] = &All.SGChemInitHCOPAbund;
      id[nt++] = REAL;

      strcpy(tag[nt], "SGChemInitHePAbund");
      addr[nt] = &All.SGChemInitHePAbund;
      id[nt++] = REAL;

      strcpy(tag[nt], "SGChemInitMPAbund");
      addr[nt] = &All.SGChemInitMPAbund;
      id[nt++] = REAL;

#ifndef SGCHEM_VARIABLE_Z
      strcpy(tag[nt], "CarbAbund");
      addr[nt] = &All.CarbAbund;
      id[nt++] = REAL;

      strcpy(tag[nt], "OxyAbund");
      addr[nt] = &All.OxyAbund;
      id[nt++] = REAL;

      strcpy(tag[nt], "MAbund");
      addr[nt] = &All.MAbund;
      id[nt++] = REAL;

      strcpy(tag[nt], "ZAtom");
      addr[nt] = &All.ZAtom;
      id[nt++] = REAL;
#endif

      strcpy(tag[nt], "DeutAbund");
      addr[nt] = &All.DeutAbund;
      id[nt++] = REAL;

      strcpy(tag[nt], "InitDustTemp");
      addr[nt] = &All.InitDustTemp;
      id[nt++] = REAL;

      strcpy(tag[nt], "UVFieldStrength");
      addr[nt] = &All.UVFieldStrength;
      id[nt++] = REAL;

      strcpy(tag[nt], "LWBGType");
      addr[nt] = &All.LWBGType;
      id[nt++] = INT;

      strcpy(tag[nt], "LWBGStartRedsh");
      addr[nt] = &All.LWBGStartRedsh;
      id[nt++] = REAL;

#ifndef SGCHEM_VARIABLE_Z
      strcpy(tag[nt], "DustToGasRatio");
      addr[nt] = &All.DustToGasRatio;
      id[nt++] = REAL;
#endif
#ifndef SGCHEM_VARIABLE_CRION
      strcpy(tag[nt], "CosmicRayIonRate");
      addr[nt] = &All.CosmicRayIonRate;
      id[nt++] = REAL;
#endif

      strcpy(tag[nt], "InitRedshift");
      addr[nt] = &All.InitRedshift;
      id[nt++] = REAL;

      strcpy(tag[nt], "ExternalDustExtinction");
      addr[nt] = &All.ExternalDustExtinction;
      id[nt++] = REAL;

      strcpy(tag[nt], "H2FormEx");
      addr[nt] = &All.H2FormEx;
      id[nt++] = REAL;

      strcpy(tag[nt], "H2FormKin");
      addr[nt] = &All.H2FormKin;
      id[nt++] = REAL;

      strcpy(tag[nt], "PhotoApprox");
      addr[nt] = &All.PhotoApprox;
      id[nt++] = INT;

      strcpy(tag[nt], "ISRFOption");
      addr[nt] = &All.ISRFOption;
      id[nt++] = INT;

      strcpy(tag[nt], "AtomicCoolOption");
      addr[nt] = &All.AtomicCoolOption;
      id[nt++] = INT;

      strcpy(tag[nt], "H2OpacityOption");
      addr[nt] = &All.H2OpacityOption;
      id[nt++] = INT;

#ifdef SGCHEM_ACCRETION_LUMINOSITY
      strcpy(tag[nt], "SGChemAccretionLuminosityOn");
      addr[nt] = &All.SGChemAccretionLuminosityOn;
      id[nt++] = INT;

      strcpy(tag[nt], "SinkAccretionRateSmoothingMass");
      addr[nt] = &All.SinkAccretionRateSmoothingMass;
      id[nt++] = REAL;
#endif

#ifdef SGCHEM_TEMPERATURE_FLOOR
      strcpy(tag[nt], "SGChemTemperatureFloor");
      addr[nt] = &All.SGChemTemperatureFloor;
      id[nt++] = REAL;
#endif

#endif /*SGCHEM*/

#ifdef TREECOLV2
      strcpy(tag[nt], "TreeColMaxDistance");
      addr[nt] = &All.TreeColMaxDistance;
      id[nt++] = REAL;
#endif

#ifdef TREECOLV2_VEL
      strcpy(tag[nt], "FracOverlap"); /* mean line overlap following Hartwig et al. 2015, ApJ, 799, 114 */
      addr[nt] = &All.FracOverlap;
      id[nt++] = REAL;
#endif

#ifdef SINKS
      strcpy(tag[nt], "NHThresh");
      addr[nt] = &SKD.NHThresh;
      id[nt++] = REAL;

      strcpy(tag[nt], "AccRad");
      addr[nt] = &SKD.AccRad;
      id[nt++] = REAL;
#endif

#ifdef SNE_FEEDBACK
      strcpy(tag[nt], "SNEInjectionCriterion");
      addr[nt] = &All.SNEInjectionCriterion;
      id[nt++] = INT;
      strcpy(tag[nt], "SNESeed");
      addr[nt] = &All.SNESeed;
      id[nt++] = INT;
      strcpy(tag[nt], "SNETargetMass");
      addr[nt] = &All.SNETargetMass;
      id[nt++] = REAL;
      strcpy(tag[nt], "SNEPeriodInYears");
      addr[nt] = &All.SNEPeriodInYears;
      id[nt++] = REAL;
      strcpy(tag[nt], "SNEMinimalParticleNumber");
      addr[nt] = &All.SNEMinimalParticleNumber;
      id[nt++] = INT;

#ifdef CLUSTERED_SNE
      strcpy(tag[nt], "SNEClusterTfin");
      addr[nt] = &All.SNEClusterTfin;
      id[nt++] = REAL;
      strcpy(tag[nt], "SNEClusterTinit");
      addr[nt] = &All.SNEClusterTinit;
      id[nt++] = REAL;
      strcpy(tag[nt], "SNENumber");
      addr[nt] = &All.SNENumber;
      id[nt++] = INT;
#endif

#ifdef INJECT_TRACER_INTO_SN
      strcpy(tag[nt], "SNETracersForEachSn");
      addr[nt] = &All.SNETracersForEachSn;
      id[nt++] = INT;
      strcpy(tag[nt], "SNETracerBitmask");
      addr[nt] = &All.SNETracerBitmask;
      id[nt++] = INT;
#endif

#if defined(SINK_PARTICLES) && defined(SINK_PARTICLES_FEEDBACK)
      strcpy(tag[nt], "SNEScatterAroundSink");
      addr[nt] = &All.SNEScatterAroundSink;
      id[nt++] = REAL;
#endif

#endif /*SNE_FEEDBACK*/

#ifdef REFINE_ONLY_WITH_TRACER
      strcpy(tag[nt], "MaxTracerVolume");
      addr[nt] = &All.MaxTracerVolume;
      id[nt++] = REAL;
      strcpy(tag[nt], "MinTracerVolume");
      addr[nt] = &All.MinTracerVolume;
      id[nt++] = REAL;
#endif

#ifdef ROTATING_HIGHRES_REGION
      strcpy(tag[nt], "Highres_x0");
      addr[nt] = &All.Highres_x0;
      id[nt++] = REAL;

      strcpy(tag[nt], "Highres_y0");
      addr[nt] = &All.Highres_y0;
      id[nt++] = REAL;

      strcpy(tag[nt], "Highres_deltaR");
      addr[nt] = &All.Highres_deltaR;
      id[nt++] = REAL;

      strcpy(tag[nt], "Highres_deltaThetaxR");
      addr[nt] = &All.Highres_deltaThetaxR;
      id[nt++] = REAL;

      strcpy(tag[nt], "Highres_deltaz");
      addr[nt] = &All.Highres_deltaz;
      id[nt++] = REAL;

      strcpy(tag[nt], "Highres_vrot");
      addr[nt] = &All.Highres_vrot;
      id[nt++] = REAL;

      strcpy(tag[nt], "Highres_t0");
      addr[nt] = &All.Highres_t0;
      id[nt++] = REAL;

      strcpy(tag[nt], "Highres_targetmass");
      addr[nt] = &All.Highres_targetmass;
      id[nt++] = REAL;
#endif

#ifdef SINK_PARTICLES
      strcpy(tag[nt], "SinkCreationDensityCodeUnits");
      addr[nt] = &SinkCreationDensityCodeUnits;
      id[nt++] = REAL;

      strcpy(tag[nt], "SinkFormationRadius");
      addr[nt] = &SinkFormationRadius;
      id[nt++] = REAL;

      strcpy(tag[nt], "SinkEvolutionDumpRateYears");
      addr[nt] = &SinkEvolutionDumpRateYears;
      id[nt++] = REAL;

#ifdef SINK_PARTICLES_FEEDBACK
      strcpy(tag[nt], "SINKStarFormationEfficiency");
      addr[nt] = &All.SINKStarFormationEfficiency;
      id[nt++] = REAL;
#endif

#ifdef SINK_PHOTOION_FEEDBACK
      strcpy(tag[nt], "StromgrenTemperature");
      addr[nt] = &All.StromgrenTemperature;
      id[nt++] = REAL;
#endif

#endif /*SINK_PARTICLES*/

#ifdef VISCOSITY
#ifdef GLOBAL_VISCOSITY
      strcpy(tag[nt], "DynamicViscosity");
      addr[nt] = &All.dyn_visc;
      id[nt++] = REAL;

      strcpy(tag[nt], "BulkViscosity");
      addr[nt] = &All.bulk_visc;
      id[nt++] = REAL;
#else
#ifdef USE_KINEMATIC_VISCOSITY
      strcpy(tag[nt], "KinematicViscosity");
      addr[nt] = &All.KinematicViscosity;
      id[nt++] = REAL;
#else
#ifdef ALPHA_VISCOSITY
      strcpy(tag[nt], "AlphaCoefficient");
      addr[nt] = &All.AlphaCoefficient;
      id[nt++] = REAL;
#endif
#endif
#endif
#endif

#ifdef THERMAL_CONDUCTION
      strcpy(tag[nt], "ThermalConductivity");
      addr[nt] = &All.ThermalConductivity;
      id[nt++] = REAL;
#endif

#ifdef TRACER_DIFFUSION
      strcpy(tag[nt], "TracerDiffusivity");
      addr[nt] = &All.TracerDiffusivity;
      id[nt++] = REAL;
#endif

#ifdef TRACER_MC
      strcpy(tag[nt], "TracerMCPerCell");
      addr[nt] = &All.TracerMCPerCell;
      id[nt++] = INT;
#endif

#ifdef TRACER_PARTICLE
      strcpy(tag[nt], "MinimumTracerHsml");
      addr[nt] = &All.MinimumTracerHsml;
      id[nt++] = REAL;
#endif

#ifdef GRAVITY_TABLE
      strcpy(tag[nt], "ExternalGravForcesFile");
      addr[nt] = All.ExternalGravForcesFile;
      id[nt++] = STRING;
#endif

#ifdef SPECIAL_BOUNDARY
      strcpy(tag[nt], "BoundaryLayerScaleFactor");
      addr[nt] = &All.BoundaryLayerScaleFactor;
      id[nt++] = REAL;

      strcpy(tag[nt], "SpecialBoundarySpeed");
      addr[nt] = &All.SpecialBoundarySpeed;
      id[nt++] = REAL;

      strcpy(tag[nt], "SpecialBoundaryMotion");
      addr[nt] = &All.SpecialBoundaryMotion;
      id[nt++] = INT;

      strcpy(tag[nt], "SpecialBoundaryType");
      addr[nt] = &All.SpecialBoundaryType;
      id[nt++] = INT;

      strcpy(tag[nt], "OutflowPressure");
      addr[nt] = &All.OutflowPressure;
      id[nt++] = REAL;
#ifdef THERMAL_CONDUCTION
      strcpy(tag[nt], "BoundaryTemperature");
      addr[nt] = &All.BoundaryTemperature;
      id[nt++] = REAL;
#endif
#endif

#if defined(COAXIAL_BOUNDARIES) && !defined(CIRCUMSTELLAR)
      strcpy(tag[nt], "InnerRadius");
      addr[nt] = &All.inner_radius;
      id[nt++] = REAL;

      strcpy(tag[nt], "OuterRadius");
      addr[nt] = &All.outer_radius;
      id[nt++] = REAL;

      strcpy(tag[nt], "InnerAngVel");
      addr[nt] = &All.omega_in;
      id[nt++] = REAL;

      strcpy(tag[nt], "OuterAngVel");
      addr[nt] = &All.omega_out;
      id[nt++] = REAL;
#endif

#ifdef AMR
      strcpy(tag[nt], "MinRefLevel");
      addr[nt] = &All.MinRefLevel;
      id[nt++] = INT;

      strcpy(tag[nt], "MaxRefLevel");
      addr[nt] = &All.MaxRefLevel;
      id[nt++] = INT;

      strcpy(tag[nt], "MeshSmoothing");
      addr[nt] = &All.AMRMeshSmoothing;
      id[nt++] = INT;
#endif

#ifdef AB_TURB
      strcpy(tag[nt], "ST_decay");
      addr[nt] = &All.StDecay;
      id[nt++] = REAL;

      strcpy(tag[nt], "ST_energy");
      addr[nt] = &All.StEnergy;
      id[nt++] = REAL;

      strcpy(tag[nt], "ST_DtFreq");
      addr[nt] = &All.StDtFreq;
      id[nt++] = REAL;

      strcpy(tag[nt], "ST_Kmin");
      addr[nt] = &All.StKmin;
      id[nt++] = REAL;

      strcpy(tag[nt], "ST_Kmax");
      addr[nt] = &All.StKmax;
      id[nt++] = REAL;

      strcpy(tag[nt], "ST_SolWeight");
      addr[nt] = &All.StSolWeight;
      id[nt++] = REAL;

      strcpy(tag[nt], "ST_AmplFac");
      addr[nt] = &All.StAmplFac;
      id[nt++] = REAL;

      strcpy(tag[nt], "ST_Seed");
      addr[nt] = &All.StSeed;
      id[nt++] = REAL;

      strcpy(tag[nt], "ST_SpectForm");
      addr[nt] = &All.StSpectForm;
      id[nt++] = INT;
#endif

#ifdef PREHEATING
      strcpy(tag[nt], "TimePreheating");
      addr[nt] = &All.TimePreheating;
      id[nt++] = REAL;

      strcpy(tag[nt], "TempPreheating");
      addr[nt] = &All.TempPreheating;
      id[nt++] = REAL;
#endif

#ifdef ADJ_BOX_POWERSPEC
      strcpy(tag[nt], "BoxWidth");
      addr[nt] = &All.BoxWidth;
      id[nt++] = REAL;

      strcpy(tag[nt], "BoxCenter_x");
      addr[nt] = &All.BoxCenter_x;
      id[nt++] = REAL;

      strcpy(tag[nt], "BoxCenter_y");
      addr[nt] = &All.BoxCenter_y;
      id[nt++] = REAL;

      strcpy(tag[nt], "BoxCenter_z");
      addr[nt] = &All.BoxCenter_z;
      id[nt++] = REAL;

      strcpy(tag[nt], "TransformSize");
      addr[nt] = &All.FourierGrid;
      id[nt++] = INT;
#endif

#ifdef VEL_POWERSPEC_BOX
      strcpy(tag[nt], "PSCenterX");
      addr[nt] = &All.PSCenterX;
      id[nt++] = REAL;

      strcpy(tag[nt], "PSCenterY");
      addr[nt] = &All.PSCenterY;
      id[nt++] = REAL;

      strcpy(tag[nt], "PSCenterZ");
      addr[nt] = &All.PSCenterZ;
      id[nt++] = REAL;

      strcpy(tag[nt], "PSRadius");
      addr[nt] = &All.PSRadius;
      id[nt++] = REAL;
#endif

#ifdef ATOMIC_DM
      strcpy(tag[nt], "ADMProtonMassInkeV");
      addr[nt] = &All.ADMProtonMassInkeV;
      id[nt++] = REAL;

      strcpy(tag[nt], "ADMElectronMassInkeV");
      addr[nt] = &All.ADMElectronMassInkeV;
      id[nt++] = REAL;

      strcpy(tag[nt], "ADMFineStructure");
      addr[nt] = &All.ADMFineStructureConstant;
      id[nt++] = REAL;
#endif

#ifdef SIDM
      strcpy(tag[nt], "DtimeFac");
      addr[nt] = &All.DtimeFac;
      id[nt++] = REAL;

      strcpy(tag[nt], "DtimeFacLim");
      addr[nt] = &All.DtimeFacLim;
      id[nt++] = REAL;

      strcpy(tag[nt], "SIDMDesNumNgb");
      addr[nt] = &All.SIDMDesNumNgb;
      id[nt++] = REAL;

      strcpy(tag[nt], "SIDMMaxNumNgbDeviation");
      addr[nt] = &All.SIDMMaxNumNgbDeviation;
      id[nt++] = REAL;
#ifdef SIDM_CONST_CROSS
      strcpy(tag[nt], "CrossSectionPerMass_in_cgs");
      addr[nt] = &All.CrossSectionPerMass_in_cgs;
      id[nt++] = REAL;
#endif
#endif

#ifdef DUST_LIVE
      strcpy(tag[nt], "DesNumNgbDust");
      addr[nt] = &All.DesNumNgbDust;
      id[nt++] = INT;

      strcpy(tag[nt], "MaxNumNgbDeviationDust");
      addr[nt] = &All.MaxNumNgbDeviationDust;
      id[nt++] = REAL;

#if !defined(DL_DRAG_SEMI_IMPLICIT) && !defined(DL_NODRAG)
      strcpy(tag[nt], "StoppingTimeFrac");
      addr[nt] = &All.StoppingTimeFrac;
      id[nt++] = REAL;
#endif

#ifdef DL_GRAIN_BINS
      strcpy(tag[nt], "MinGrainSize");
      addr[nt] = &All.MinGrainSize;
      id[nt++] = REAL;

      strcpy(tag[nt], "MaxGrainSize");
      addr[nt] = &All.MaxGrainSize;
      id[nt++] = REAL;

      strcpy(tag[nt], "GrainDensity");
      addr[nt] = &All.GrainDensity;
      id[nt++] = REAL;

      strcpy(tag[nt], "MaxBinFracMassChg");
      addr[nt] = &All.MaxBinFracMassChg;
      id[nt++] = REAL;

#if defined(DL_SHATTERING) || defined(DL_COAGULATION) || defined(DL_PRODUCTION)
      strcpy(tag[nt], "GrainDataPath");
      addr[nt] = &All.GrainDataPath;
      id[nt++] = STRING;
#endif

#ifdef DL_PRODUCTION
      strcpy(tag[nt], "DustTargetFrac");
      addr[nt] = &All.DustTargetFrac;
      id[nt++] = REAL;

      strcpy(tag[nt], "NumDustPerSpawn");
      addr[nt] = &All.NumDustPerSpawn;
      id[nt++] = INT;

#ifdef DL_REFINEMENT
      strcpy(tag[nt], "DustMaxFrac");
      addr[nt] = &All.DustMaxFrac;
      id[nt++] = REAL;
#endif

#ifdef DL_DEREFINEMENT
      strcpy(tag[nt], "DustMinFrac");
      addr[nt] = &All.DustMinFrac;
      id[nt++] = REAL;
#endif
#endif /* DL_PRODUCTION */

#ifdef DL_SUBCYCLE
      strcpy(tag[nt], "DustSubcycleFac");
      addr[nt] = &All.DustSubcycleFac;
      id[nt++] = REAL;
#endif
#endif /* DL_GRAIN_BINS */

#ifdef DL_RADIATION
      strcpy(tag[nt], "GrainRPPath");
      addr[nt] = &All.GrainRPPath;
      id[nt++] = STRING;
#endif
#endif /* DUST_LIVE */

#ifdef REFINEMENT_VOLUME_LIMIT
      strcpy(tag[nt], "MaxVolumeDiff");
      addr[nt] = &All.MaxVolumeDiff;
      id[nt++] = REAL;

      strcpy(tag[nt], "MinVolume");
      addr[nt] = &All.MinVolume;
      id[nt++] = REAL;

      strcpy(tag[nt], "MaxVolume");
      addr[nt] = &All.MaxVolume;
      id[nt++] = REAL;
#endif

#ifdef REFINEMENT_BY_DENSITY
      strcpy(tag[nt], "MinimumDensityForRefinement");
      addr[nt] = &All.MinimumDensityForRefinement;
      id[nt++] = REAL;

      strcpy(tag[nt], "MinimumVolumeForDensityRefinement");
      addr[nt] = &All.MinimumVolumeForDensityRefinement;
      id[nt++] = REAL;
#endif

#ifdef TILE_ICS
      strcpy(tag[nt], "TileICsFactor");
      addr[nt] = &All.TileICsFactor;
      id[nt++] = INT;
#endif

#ifdef OTVET
      strcpy(tag[nt], "IonizingLumPerSolarMass");
      addr[nt] = &All.IonizingLumPerSolarMass;
      id[nt++] = REAL;
#ifdef OTVET_MULTI_FREQUENCY
      strcpy(tag[nt], "star_Teff");
      addr[nt] = &All.star_Teff;
      id[nt++] = REAL;
#endif
#ifdef OTVET_SCATTER_SOURCE
      strcpy(tag[nt], "OtvetNumNgbSource");
      addr[nt] = &All.OtvetNumNgbSource;
      id[nt++] = REAL;
      strcpy(tag[nt], "OtvetMaxNumNgbDeviationSource");
      addr[nt] = &All.OtvetMaxNumNgbDeviationSource;
      id[nt++] = REAL;
#endif
#endif

#if defined(CONDUCTION) || defined(MONOTONE_CONDUCTION)
      strcpy(tag[nt], "MaxConductionSubCycles");
      addr[nt] = &All.MaxConductionSubCycles;
      id[nt++] = INT;

      strcpy(tag[nt], "ConductionEfficiency");
      addr[nt] = &All.ConductionEfficiency;
      id[nt++] = REAL;

#ifdef RESTRICT_KAPPA
      strcpy(tag[nt], "MaxDiffusivity");
      addr[nt] = &All.MaxDiffusivity;
      id[nt++] = REAL;
#endif

#endif

#ifdef RADPRESS_OPT_THIN
      /*
            strcpy(tag[nt], "IonizingLumPerSolarMass");
            addr[nt] = &All.IonizingLumPerSolarMass;
            id[nt++] = REAL;
      */
      strcpy(tag[nt], "RadPressure_MaxAge");
      addr[nt] = &All.RadPressure_MaxAge;
      id[nt++] = REAL;
#endif

#ifdef ADDBACKGROUNDGRID
      strcpy(tag[nt], "GridSize");
      addr[nt] = &All.GridSize;
      id[nt++] = INT;
#endif

#if defined(SHOCK_FINDER_BEFORE_OUTPUT) || defined(SHOCK_FINDER_ON_THE_FLY)
      // standard parameters
      All.RayMemFac   = 1.5;
      All.MachMin     = 1.3;
      All.RayStepsMax = 5;
#endif

#ifdef SHOCK_FINDER_POST_PROCESSING
      strcpy(tag[nt], "RayMemFac");
      addr[nt] = &All.RayMemFac;
      id[nt++] = REAL;

      strcpy(tag[nt], "MachMin");
      addr[nt] = &All.MachMin;
      id[nt++] = REAL;

      strcpy(tag[nt], "RayStepsMax");
      addr[nt] = &All.RayStepsMax;
      id[nt++] = INT;

      strcpy(tag[nt], "NumFilesPerOutput");
      addr[nt] = &All.NumFilesPerOutput;
      id[nt++] = INT;

      strcpy(tag[nt], "OutputDirShockFinder");
      addr[nt] = All.OutputDirShockFinder;
      id[nt++] = STRING;
#endif

#if defined(SHOCK_FINDER_POST_PROCESSING) || defined(SHOCK_FINDER_BEFORE_OUTPUT) || defined(SHOCK_FINDER_ON_THE_FLY)
#ifdef SKIP_BORDER
      strcpy(tag[nt], "DistToBorder");
      addr[nt] = &All.DistToBorder;
      id[nt++] = REAL;
#endif
#ifdef SHOCK_SUBVOLUME
      strcpy(tag[nt], "ShockSubvolumeNum");
      addr[nt] = &All.ShockSubvolumeNum;
      id[nt++] = INT;
#endif
#endif

#ifdef DG
#if defined(CHARACTERISTIC_LIMITER) || defined(CONSERVED_LIMITER)
      strcpy(tag[nt], "DG_beta");
      addr[nt] = &All.DG_beta;
      id[nt++] = REAL;
#endif

#ifdef DISCONTINUITY_DETECTION
      strcpy(tag[nt], "DG_alpha");
      addr[nt] = &All.DG_alpha;
      id[nt++] = REAL;
#endif

#ifdef ANGLE_BOUND
      strcpy(tag[nt], "DG_min_angle");
      addr[nt] = &All.DG_min_angle;
      id[nt++] = REAL;
#endif

#ifdef MINMOD_B
      strcpy(tag[nt], "DG_M");
      addr[nt] = &All.DG_M;
      id[nt++] = REAL;
#endif

#ifdef MACHNUM_B
      strcpy(tag[nt], "DG_lim_mach_min");
      addr[nt] = &All.DG_lim_mach_min;
      id[nt++] = REAL;
#endif

#if defined(REFINEMENT_SPLIT_CELLS) || defined(REFINEMENT_MERGE_CELLS)
      strcpy(tag[nt], "DG_TargetSlope");
      addr[nt] = &All.DG_TargetSlope;
      id[nt++] = REAL;

      strcpy(tag[nt], "DG_SlopeRangeFactor");
      addr[nt] = &All.DG_SlopeRangeFactor;
      id[nt++] = REAL;
#endif

#ifdef RK2
      All.DG_RK2_alpha = 1;  // use SSP RK2 Heun's method as a standard
#endif
#endif

#ifdef COSMIC_RAYS
      strcpy(tag[nt], "GammaCR");
      addr[nt] = &All.GammaCR;
      id[nt++] = REAL;

#ifdef COSMIC_RAYS_SN_INJECTION
      strcpy(tag[nt], "CREnergyInputPerSolarMassOfStarFormation");
      addr[nt] = &All.CREnergyInputPerSolarMassOfStarFormation;
      id[nt++] = REAL;

      strcpy(tag[nt], "DesNumNgbCRInjection");
      addr[nt] = &All.DesNumNgbCRInjection;
      id[nt++] = REAL;

      strcpy(tag[nt], "MaxNumNgbDeviationCRInjection");
      addr[nt] = &All.MaxNumNgbDeviationCRInjection;
      id[nt++] = REAL;
#endif

      strcpy(tag[nt], "MinimumCREnergyDensity");
      addr[nt] = &All.MinimumCREnergyDensity;
      id[nt++] = REAL;

#ifdef COSMIC_RAYS_SHOCK_ACCELERATION
      strcpy(tag[nt], "AccelerationEfficiency");
      addr[nt] = &All.AccelerationEfficiency;
      id[nt++] = REAL;

      strcpy(tag[nt], "CriticalMachnumber");
      addr[nt] = &All.CriticalMachnumber;
      id[nt++] = REAL;
#endif

#ifdef COSMIC_RAYS_DIFFUSION
#ifdef COSMIC_RAYS_DIFFUSION_CONSTANT_TIMESTEP
      strcpy(tag[nt], "CRDiffusionTimestep");
      addr[nt] = &All.CRDiffusionTimestep;
      id[nt++] = REAL;
#else
      strcpy(tag[nt], "CourantCRDiffusion");
      addr[nt] = &All.CourantCRDiffusion;
      id[nt++] = REAL;
#endif
      strcpy(tag[nt], "CRDiffusionCoefficient");
      addr[nt] = &All.CR_Diffusion_Coefficient;
      id[nt++] = REAL;

#ifdef COSMIC_RAYS_DIFFUSION_BOUNDARY_X
      strcpy(tag[nt], "CR_Boundary_X_Upper");
      addr[nt] = &All.CR_Boundary_X_Upper;
      id[nt++] = REAL;

      strcpy(tag[nt], "CR_Boundary_X_Lower");
      addr[nt] = &All.CR_Boundary_X_Lower;
      id[nt++] = REAL;
#endif

#ifdef COSMIC_RAYS_DIFFUSION_BOUNDARY_Y
      strcpy(tag[nt], "CR_Boundary_Y_Upper");
      addr[nt] = &All.CR_Boundary_Y_Upper;
      id[nt++] = REAL;

      strcpy(tag[nt], "CR_Boundary_Y_Lower");
      addr[nt] = &All.CR_Boundary_Y_Lower;
      id[nt++] = REAL;
#endif

#ifdef COSMIC_RAYS_DIFFUSION_BOUNDARY_Z
      strcpy(tag[nt], "CR_Boundary_Z_Upper");
      addr[nt] = &All.CR_Boundary_Z_Upper;
      id[nt++] = REAL;

      strcpy(tag[nt], "CR_Boundary_Z_Lower");
      addr[nt] = &All.CR_Boundary_Z_Lower;
      id[nt++] = REAL;
#endif

#endif

#ifdef COSMIC_RAYS_STREAMING_EXPLICIT
      strcpy(tag[nt], "CR_Chi");
      addr[nt] = &All.CR_Chi;
      id[nt++] = REAL;
#endif

#endif

#ifdef FLD
#ifdef FLD_TEST_BOUNDARY
      strcpy(tag[nt], "FLDFlux");
      addr[nt] = &All.fld_Flux;
      id[nt++] = REAL;
      strcpy(tag[nt], "FLDGamma");
      addr[nt] = &All.fld_n_gamma;
      id[nt++] = REAL;
      strcpy(tag[nt], "FLDDensity");
      addr[nt] = &All.fld_density;
      id[nt++] = REAL;
      strcpy(tag[nt], "FLDEgySpec");
      addr[nt] = &All.fld_u;
      id[nt++] = REAL;
#endif
#ifdef FLD_CONST_KAPPA
      strcpy(tag[nt], "Kappa_R");
      addr[nt] = &All.Kappa_R;
      id[nt++] = REAL;
      strcpy(tag[nt], "Kappa_P");
      addr[nt] = &All.Kappa_P;
      id[nt++] = REAL;
#endif
#ifdef FLD_MARSHAK
      strcpy(tag[nt], "Epsilon");
      addr[nt] = &All.Epsilon;
      id[nt++] = REAL;
#endif
#endif

#ifdef EXTERNALSHEARBOX
      strcpy(tag[nt], "ShearBoxSigma0");
      addr[nt] = &All.ShearBoxSigma0;
      id[nt++] = REAL;

      strcpy(tag[nt], "ShearBoxFg");
      addr[nt] = &All.ShearBoxFg;
      id[nt++] = REAL;

      strcpy(tag[nt], "ShearBoxMu");
      addr[nt] = &All.ShearBoxMu;
      id[nt++] = REAL;

#endif

#ifdef LOCAL_FEEDBACK
      strcpy(tag[nt], "LocalFeedbackSNEnergy");
      addr[nt] = &All.LocalFeedbackSNEnergy;
      id[nt++] = REAL;

      strcpy(tag[nt], "LocalFeedbackSNMassReturn");
      addr[nt] = &All.LocalFeedbackSNMassReturn;
      id[nt++] = REAL;

      strcpy(tag[nt], "LocalFeedbackSNRate");
      addr[nt] = &All.LocalFeedbackSNRate;
      id[nt++] = REAL;

#if defined(COSMIC_RAYS) && !defined(COSMIC_RAYS_SHOCK_ACCELERATION)
      strcpy(tag[nt], "LocalFeedbackCRInjectionFraction");
      addr[nt] = &All.LocalFeedbackCRInjectionFraction;
      id[nt++] = REAL;
#endif

#ifdef LOCAL_FEEDBACK_PARTICLES
      strcpy(tag[nt], "LocalFeedbackSNTimeDelay");
      addr[nt] = &All.LocalFeedbackSNTimeDelay;
      id[nt++] = REAL;

      strcpy(tag[nt], "LocalFeedbackSNTimeSpread");
      addr[nt] = &All.LocalFeedbackSNTimeSpread;
      id[nt++] = REAL;

#endif

#if !defined(EXTERNALSHEARBOX_KSRATE) || defined(EXTERNALSHEARBOX_MIXED_INJECTION)
      strcpy(tag[nt], "LocalFeedbackSFEff");
      addr[nt] = &All.LocalFeedbackSFEff;
      id[nt++] = REAL;

      strcpy(tag[nt], "LocalFeedbackSFDenThresh");
      addr[nt] = &All.LocalFeedbackSFDenThresh;
      id[nt++] = REAL;

#endif

#ifdef LOCAL_KINETIC
      strcpy(tag[nt], "LocalFeedbackKineticInjectionFraction");
      addr[nt] = &All.LocalFeedbackKineticInjectionFraction;
      id[nt++] = REAL;
#endif
#endif

#ifdef LOCALIZED_SOFTENINGS
      strcpy(tag[nt], "LocalizedSofteningsParticleNumber");
      addr[nt] = &All.LocalizedSofteningsParticleNumber;
      id[nt++] = INT;
#endif

#ifdef NON_IDEAL_MHD
#if defined(OHMIC_DIFFUSION) || defined(IMPLICIT_OHMIC_DIFFUSION)
      strcpy(tag[nt], "OhmicDiffusionCoefficient");
      addr[nt] = &All.OhmicDiffusionCoefficient;
      id[nt++] = REAL;
#endif
#ifdef AMBIPOLAR_DIFFUSION
      strcpy(tag[nt], "AmbipolarDiffusionCoefficient");
      addr[nt] = &All.AmbipolarDiffusionCoefficient;
      id[nt++] = REAL;
#endif
#ifdef NON_IDEAL_MHD_EXPLICIT_LIMIT_TIMESTEP
      strcpy(tag[nt], "NonidealMHDTimelimit");
      addr[nt] = &All.NonidealMHDTimelimit;
      id[nt++] = REAL;
#endif
#endif

#ifdef SF_STELLAR_MASS_TO_GAS_MASS_RATIO
      strcpy(tag[nt], "StellarMassToGasMassRatio");
      addr[nt] = &All.StellarMassToGasMassRatio;
      id[nt++] = REAL;
#endif

#ifdef AURIGA_MOVIE
      strcpy(tag[nt], "Auriga_Movie_CenterRadius");
      addr[nt] = &All.Auriga_Movie_CenterRadius;
      id[nt++] = REAL;

      strcpy(tag[nt], "Auriga_Movie_Directory");
      addr[nt] = All.Auriga_Movie_Directory;
      id[nt++] = STRING;

      int ix, iy;
      for(ix = 0; ix < 3; ix++)
        for(iy = 0; iy < 3; iy++)
          {
            sprintf(buf, "Auriga_Movie_Galaxy_Rotation%d%d", ix, iy);
            strcpy(tag[nt], buf);
            addr[nt] = &All.Auriga_Movie_Galaxy_Rotation[ix][iy];
            id[nt++] = REAL;
          }

      strcpy(tag[nt], "Auriga_Movie_OutputListFilename");
      addr[nt] = All.Auriga_Movie_OutputListFilename;
      id[nt++] = STRING;
#endif

#ifdef HCOUTPUT
      strcpy(tag[nt], "HCOutput_RadialCut");
      addr[nt] = &All.HCOutput_RadialCut;
      id[nt++] = REAL;

      strcpy(tag[nt], "HCOutput_CenterRadius");
      addr[nt] = &All.HCOutput_CenterRadius;
      id[nt++] = REAL;

      strcpy(tag[nt], "HCOutput_Directory");
      addr[nt] = All.HCOutput_Directory;
      id[nt++] = STRING;

      strcpy(tag[nt], "HCOutput_OutputListFilename");
      addr[nt] = All.HCOutput_OutputListFilename;
      id[nt++] = STRING;
#endif

#ifdef ONEDIMS_SPHERICAL
      strcpy(tag[nt], "CoreRadius");
      addr[nt] = &All.CoreRadius;
      id[nt++] = REAL;

      strcpy(tag[nt], "CoreMass");
      addr[nt] = &All.CoreMass;
      id[nt++] = REAL;
#endif

#ifdef MODGRAV
      strcpy(tag[nt], "ModifiedGravityFile");
      addr[nt] = &All.ModifiedGravityFile;
      id[nt++] = REAL;
      strcpy(tag[nt], "MaxAMRLevel");
      addr[nt] = &All.MaxAMRLevel;
      id[nt++] = REAL;
      strcpy(tag[nt], "MinLevelTopLeaf");
      addr[nt] = &All.MinLevelTopLeaf;
      id[nt++] = REAL;
      strcpy(tag[nt], "MaxLevelTopLeaf");
      addr[nt] = &All.MaxLevelTopLeaf;
      id[nt++] = REAL;
      strcpy(tag[nt], "MaxLevelFullTree");
      addr[nt] = &All.MaxLevelFullTree;
      id[nt++] = REAL;
#endif

#ifdef SGS_TURBULENCE
      strcpy(tag[nt], "MinimumSgsTSpecificEnergy");
      addr[nt] = &All.SgsTConst.MinimumSgsTSpecificEnergy;
      id[nt++] = REAL;

#ifdef SGS_TURBULENCE_EDDY_VISCOSITY
      strcpy(tag[nt], "EddyViscosityParameterCnu");
      addr[nt] = &All.SgsTConst.EddyViscosityParameterCnu;
      id[nt++] = REAL;

      strcpy(tag[nt], "EddyViscosityParameterCepsilon");
      addr[nt] = &All.SgsTConst.EddyViscosityParameterCepsilon;
      id[nt++] = REAL;
#endif /* #ifdef SGS_TURBULENCE_EDDY_VISCOSITY */
#endif /* #ifdef SGS_TURBULENCE */

#ifdef SIMPLEX
#if SX_SOURCES == 10
      strcpy(tag[nt], "TestSrcFile");
      addr[nt] = All.sxTestSrcFile;
      id[nt++] = STRING;
#endif
      strcpy(tag[nt], "UnitPhotons_per_s");
      addr[nt] = &All.UnitPhotons_per_s;
      id[nt++] = REAL;
      strcpy(tag[nt], "MinNumPhotons");
      addr[nt] = &All.sxMinNumPhotons;
      id[nt++] = REAL;
#endif

#ifdef AGB_WIND
      strcpy(tag[nt], "AGBWindDensity");
      addr[nt] = &All.AGBWindDensity;
      id[nt++] = REAL;

      strcpy(tag[nt], "AGBWindVelocity");
      addr[nt] = &All.AGBWindVelocity;
      id[nt++] = REAL;

      strcpy(tag[nt], "AGBWindSpecificEnergy");
      addr[nt] = &All.AGBWindSpecificEnergy;
      id[nt++] = REAL;

      strcpy(tag[nt], "AGBWindCenterX");
      addr[nt] = &All.AGBWindCenterX;
      id[nt++] = REAL;

      strcpy(tag[nt], "AGBWindCenterY");
      addr[nt] = &All.AGBWindCenterY;
      id[nt++] = REAL;

      strcpy(tag[nt], "AGBWindCenterZ");
      addr[nt] = &All.AGBWindCenterZ;
      id[nt++] = REAL;
#endif

      if((fd = fopen(fname, "r")))
        {
          sprintf(buf, "%s%s", fname, "-usedvalues");
          if(!(fdout = fopen(buf, "w")))
            {
              printf("error opening file '%s' \n", buf);
              errorFlag = 1;
            }
          else
            {
              printf("Obtaining parameters from file '%s':\n\n", fname);
              while(!feof(fd))
                {
                  *buf         = 0;
                  char *status = fgets(buf, MAXLEN_PARAM_TAG + MAXLEN_PARAM_VALUE + 200, fd);
                  if(sscanf(buf, "%s%s%s", buf1, buf2, buf3) < 2)
                    continue;

                  if(buf1[0] == '%')
                    continue;

                  for(i = 0, j = -1; i < nt; i++)
                    if(strcmp(buf1, tag[i]) == 0)
                      {
                        if(param_handled[i] == 0)
                          {
                            j                = i;
                            param_handled[i] = 1;
                            break;
                          }
                        else
                          {
                            j = -2;
                            break;
                          }
                      }

                  if(j >= 0)
                    {
                      switch(id[j])
                        {
                          case REAL:
                            *((double *)addr[j]) = atof(buf2);
                            sprintf(buf3, "%%-%ds%%g\n", MAXLEN_PARAM_TAG);
                            fprintf(fdout, buf3, buf1, *((double *)addr[j]));
                            fprintf(stdout, "        ");
                            fprintf(stdout, buf3, buf1, *((double *)addr[j]));
                            break;
                          case STRING:
                            strcpy((char *)addr[j], buf2);
                            sprintf(buf3, "%%-%ds%%s\n", MAXLEN_PARAM_TAG);
                            fprintf(fdout, buf3, buf1, buf2);
                            fprintf(stdout, "        ");
                            fprintf(stdout, buf3, buf1, buf2);
                            break;
                          case INT:
                            *((int *)addr[j]) = atoi(buf2);
                            sprintf(buf3, "%%-%ds%%d\n", MAXLEN_PARAM_TAG);
                            fprintf(fdout, buf3, buf1, *((int *)addr[j]));
                            fprintf(stdout, "        ");
                            fprintf(stdout, buf3, buf1, *((int *)addr[j]));
                            break;
                        }
                    }
                  else if(j == -2)
                    {
#ifdef ALLOWEXTRAPARAMS
                      warn("Tag '%s' ignored from file %s !", buf1, fname);
#else
                      fprintf(stdout, "Error in file %s:   Tag '%s' multiply defined.\n", fname, buf1);
                      errorFlag = 1;
#endif
                    }
                  else
                    {
#ifdef ALLOWEXTRAPARAMS
                      warn("Tag '%s' ignored from file %s !", buf1, fname);
#else
                      fprintf(stdout, "Error in file %s:   Tag '%s' not allowed\n", fname, buf1);
                      errorFlag = 1;
#endif
                    }
                }
              fclose(fd);
              fclose(fdout);
              printf("\n");

              i = strlen(All.OutputDir);
              if(i > 0)
                if(All.OutputDir[i - 1] != '/')
                  strcat(All.OutputDir, "/");

#ifdef SHOCK_FINDER_POST_PROCESSING
              mkdir(All.OutputDirShockFinder, 02755);
              sprintf(buf1, "%s%s", fname, "-usedvalues");
              sprintf(buf2, "%s/%s", All.OutputDirShockFinder, "parameters-usedvalues");
              printf("%s\n", buf2);
#ifndef NOCALLSOFSYSTEM
              if(errorFlag == 0)
                system(buf3);
#endif
#else
              mkdir(All.OutputDir, 02755);
              sprintf(buf1, "%s%s", fname, "-usedvalues");
              sprintf(buf2, "%s%s", All.OutputDir, "parameters-usedvalues");
              sprintf(buf3, "cp %s %s", buf1, buf2);
#ifndef NOCALLSOFSYSTEM
              if(errorFlag == 0)
                {
                  int status = system(buf3);
                }
#endif
#endif
            }
        }
      else
        {
          printf("Parameter file %s not found.\n", fname);
          errorFlag = 1;
        }

      for(i = 0; i < nt; i++)
        {
          if(param_handled[i] != 1)
            {
              printf("Error. I miss a value for tag '%s' in parameter file '%s'.\n", tag[i], fname);
              errorFlag = 1;
            }
        }

      if(All.OutputListOn && errorFlag == 0)
        errorFlag += read_outputlist(All.OutputListFilename);
      else
        All.OutputListLength = 0;
    }

  MPI_Bcast(&errorFlag, 1, MPI_INT, 0, MPI_COMM_WORLD);

  if(errorFlag)
    {
      MPI_Finalize();
      exit(errorFlag);
    }

  All.NParameters = nt;

  /* now communicate the relevant parameters to the other processes */
  MPI_Bcast(&All, sizeof(struct global_data_all_processes), MPI_BYTE, 0, MPI_COMM_WORLD);
#ifdef CHIMES
  MPI_Bcast(&ChimesGlobalVars, sizeof(struct globalVariables), MPI_BYTE, 0, MPI_COMM_WORLD);
#endif

#ifdef TOLERATE_WRITE_ERROR
  MPI_Bcast(AlternativeOutputDir, MAXLEN_PATH, MPI_BYTE, 0, MPI_COMM_WORLD);
#endif

#ifdef HOST_MEMORY_REPORTING
  check_maxmemsize_setting();
#endif

  mymalloc_init();

  Parameters      = (char(*)[MAXLEN_PARAM_TAG])mymalloc("Parameters", All.NParameters * MAXLEN_PARAM_TAG * sizeof(char));
  ParametersValue = (char(*)[MAXLEN_PARAM_VALUE])mymalloc("ParametersValue", All.NParameters * MAXLEN_PARAM_VALUE * sizeof(char));
  ParametersType  = (char *)mymalloc("ParamtersType", All.NParameters * sizeof(char));

  if(ThisTask == 0)
    {
      for(i = 0; i < All.NParameters; i++)
        {
          strncpy(Parameters[i], tag[i], MAXLEN_PARAM_TAG);
          ParametersType[i] = id[i];
          void *tmp         = ParametersValue[i];
          switch(id[i])
            {
              case REAL:
                *((double *)tmp) = *((double *)addr[i]);
                break;
              case STRING:
                strncpy((char *)tmp, (char *)addr[i], MAXLEN_PARAM_VALUE);
                break;
              case INT:
                tmp           = ParametersValue[i];
                *((int *)tmp) = *((int *)addr[i]);
                break;
            }
        }
    }

  MPI_Bcast(Parameters, sizeof(char) * All.NParameters * MAXLEN_PARAM_TAG, MPI_BYTE, 0, MPI_COMM_WORLD);
  MPI_Bcast(ParametersValue, sizeof(char) * All.NParameters * MAXLEN_PARAM_VALUE, MPI_BYTE, 0, MPI_COMM_WORLD);
  MPI_Bcast(ParametersType, sizeof(char) * All.NParameters, MPI_BYTE, 0, MPI_COMM_WORLD);

  for(i = 0; i < NTYPES; i++)
    {
      if(All.SofteningTypeOfPartType[i] >= NSOFTTYPES || All.SofteningTypeOfPartType[i] < 0)
        {
          mpi_printf("SofteningTypeOfPartType%  invalid (NSOFTTYPES=%d)\n", i, NSOFTTYPES);
          errorFlag = 1;
        }
    }

  if(errorFlag)
    mpi_terminate("Softening invalid!");

  if(All.NumFilesWrittenInParallel > NTask)
    {
      if(ThisTask == 0)
        warn("NOTICE: Reducing requested NumFilesWrittenInParallel=%d to %d\n", All.NumFilesWrittenInParallel, NTask);
      All.NumFilesWrittenInParallel = NTask;
    }

  if(All.NumFilesWrittenInParallel == 0)
    {
      mpi_printf("NOTICE: All.NumFilesWrittenInParallel has been set to be equal to the number of processors\n");
      All.NumFilesWrittenInParallel = NTask;
    }

#ifdef SUBBOX_SNAPSHOTS
  if(All.SubboxNumFilesWrittenInParallel > NTask)
    {
      if(ThisTask == 0)
        warn("NOTICE: Reducing requested SubboxNumFilesWrittenInParallel=%d to %d\n", All.SubboxNumFilesWrittenInParallel, NTask);
      All.SubboxNumFilesWrittenInParallel = NTask;
    }
  if(All.SubboxNumFilesWrittenInParallel > All.SubboxNumFilesPerSnapshot)
    {
      if(ThisTask == 0)
        warn("Reducing requested SubboxNumFilesWrittenInParallel=%d to be equal SubboxNumFilesPerSnapshot=%d",
             All.SubboxNumFilesWrittenInParallel, All.SubboxNumFilesPerSnapshot);
      All.SubboxNumFilesWrittenInParallel = All.SubboxNumFilesPerSnapshot;
    }
#endif

#ifndef GRAVITY_NOT_PERIODIC
  if(All.PeriodicBoundariesOn == 0)
    {
      mpi_terminate(
          "Code was compiled with gravity periodic boundary conditions switched on.\nYou must set `PeriodicBoundariesOn=1', or "
          "recompile the code.\n");
    }
#else
  if(All.PeriodicBoundariesOn == 1)
    {
      mpi_terminate(
          "Code was compiled with gravity periodic boundary conditions switched off.\nYou must set `PeriodicBoundariesOn=0', or "
          "recompile the code.\n");
    }
#endif

#ifdef COOLING
  if(All.CoolingOn == 0)
    {
      mpi_terminate("Code was compiled with cooling switched on.\nYou must set `CoolingOn=1', or recompile the code.\n");
    }
#else
  if(All.CoolingOn == 1)
    {
      mpi_terminate("Code was compiled with cooling switched off.\nYou must set `CoolingOn=0', or recompile the code.\n");
    }
#endif

  if(All.TypeOfTimestepCriterion >= 3)
    {
      mpi_terminate("The specified timestep criterion\nis not valid\n");
    }

    // CLEANUP move some of the printf statements to #error directives
#if(NTYPES < 6)
  mpi_terminate("NTYPES < 6 is not allowed.\n");
#endif

#if(NTYPES > 15)
  mpi_terminate("NTYPES > 15 is not supported yet.\n");
#endif

#if(NTYPES > 8)
  if(All.ICFormat == SNAP_FORMAT_GADGET || All.ICFormat == SNAP_FORMAT_GADGET_VARIANT)
    {
      mpi_terminate("NTYPES>8 is not allowed with ICFormat=%d, since the header block is limited to 256 bytes.\n", All.ICFormat);
    }
#endif

#ifdef USE_SFR
  if(All.StarformationOn == 0)
    {
      mpi_terminate("Code was compiled with star formation switched on.\nYou must set `StarformationOn=1', or recompile the code.\n");
    }
#ifndef LOCAL_FEEDBACK
  if(All.CoolingOn == 0)
    {
      mpi_terminate(
          "You try to use the code with star formation enabled,\nbut you did not switch on cooling.\nThis mode is not supported.\n");
    }
#endif
#else
  if(All.StarformationOn == 1)
    {
      mpi_terminate("Code was compiled with star formation switched off.\nYou must set `StarformationOn=0', or recompile the code.\n");
    }
#endif

    /*
    #if defined(VORONOI_STATIC_MESH) && defined(VORONOI_DYNAMIC_UPDATE)
      mpi_terminate("Code was compiled with VORONOI_STATIC_MESH and VORONOI_DYNAMIC_UPDATE.\n This is not allowed!\n");
    #endif
    */

#ifdef TRACER_PARTICLE
#if((TRACER_PART_TMAX_TIME) || (TRACER_PART_TMAX_RHO)) && !(TRACER_PART_TMAX)
  mpi_terminate(
      "Code was compiled with TRACER_PART_TMAX_TIME or TRACER_PART_TMAX_RHO.\n You must compile with TRACER_PART_TMAX as well.\n");
#endif
#endif

#if defined(TRACER_PART_NUM_FLUID_QUANTITIES) && !defined(TRACER_PARTICLE)
  mpi_terminate("Code was compiled with TRACER_PART_NUM_FLUID_QUANTITIES but without TRACER_PARTICLE.\n This is not allowed.\n");
#endif

#if defined(TRACER_MC_NUM_FLUID_QUANTITIES) && !defined(TRACER_MC)
  mpi_terminate("Code was compiled with TRACER_MC_NUM_FLUID_QUANTITIES but without TRACER_MC.\n This is not allowed.\n");
#endif

#ifdef TRACER_MC
#ifdef REFINEMENT_MERGE_PAIRS
  mpi_terminate("Code was compiled with TRACER_MC and REFINEMENT_MERGE_PAIRS together.\n This is not supported yet.\n");
#endif
#if((TRACER_MC_TMAX_TIME) || (TRACER_MC_TMAX_RHO)) && !(TRACER_MC_TMAX)
  mpi_terminate("Code was compiled with TRACER_MC_TMAX_TIME or TRACER_MC_TMAX_RHO.\n You must compile with TRACER_MC_TMAX as well.\n");
#endif
#if TRACER_MC_LAST_STAR_TIME && !(defined(GFM_STELLAR_EVOLUTION) || defined(GFM_WINDS))
  mpi_terminate(
      "Code was compiled with TRACER_MC and TRACER_MC_LAST_STAR_TIME but without GFM_STELLAR_EVOLUTION or GFM_WINDS.\n There is no "
      "point doing that.\n");
#endif
#if TRACER_MC_WIND_COUNTER && !defined(GFM_WINDS)
  mpi_terminate("Code was compiled with TRACER_MC_WIND_COUNTER but without GFM_WINDS.\n This is not allowed.\n");
#endif
#if defined(ADD_GROUP_PROPERTIES) || defined(RECOMPUTE_POTENTIAL_IN_SNAPSHOT)
  mpi_terminate("TRACER_MC compiled with ADD_GROUP_PROPERTIES or RECOMPUTE_POTENTIAL_IN_SNAPSHOT. More thought required in I/O.\n");
#endif
#endif /* TRACER_MC */

#if defined(ENFORCE_JEANS_STABILITY_OF_CELLS) &&                                                              \
    (defined(ISOTHERM_EQS) || defined(LOCALLY_ISOTHERM_DISK) || defined(TGCHEM) || defined(EOS_DEGENERATE) || \
     (defined(USE_SFR) && !defined(FM_SFR)) || defined(EOS_OPAL))
  if(ThisTask == 0)
    warn("Code was compiled with ENFORCE_JEANS_STABILITY_OF_CELLS together with another EOS. Please make sure you really want this.");
#endif

#if(defined(GFM_CONST_IMF) || defined(GFM_VARIABLE_IMF)) && !defined(GFM_STELLAR_EVOLUTION)
  mpi_terminate("Code was compiled with GFM_*_IMF, but not with GFM_STELLAR_EVOLUTION.\nThis is not allowed.\n");
#endif

#if defined(GFM_DUST) && !defined(GFM_STELLAR_EVOLUTION)
  mpi_terminate("Code was compiled with GFM_DUST, but not with GFM_STELLAR_EVOLUTION.\nThis is not allowed.\n");
#endif

#if defined(DL_GRAIN_BINS) && !defined(GFM_STELLAR_EVOLUTION)
  mpi_terminate("Code was compiled with DL_GRAIN_BINS, but not with GFM_STELLAR_EVOLUTION.\nThis is not allowed.\n");
#endif

#if defined(GFM_CHEMTAGS) && !defined(GFM_STELLAR_EVOLUTION)
  mpi_terminate("Code was compiled with GFM_CHEMTAGS, but not with GFM_STELLAR_EVOLUTION.\nThis is not allowed.\n");
#endif

#if defined(GFM_SPLITFE) && !defined(GFM_STELLAR_EVOLUTION)
  mpi_terminate("Code was compiled with GFM_SPLITFE, but not with GFM_STELLAR_EVOLUTION.\nThis is not allowed.\n");
#endif

#if defined(GFM_SPLITFE) && !defined(GFM_CHEMTAGS)
  mpi_terminate("Code was compiled with GFM_SPLITFE, but not with GFM_CHEMTAGS.\nThis is not allowed.\n");
#endif

#if defined(GFM_RPROCESS) && !defined(GFM_SPLITFE)
  mpi_terminate("Code was compiled with GFM_RPROCESS, but not with GFM_SPLITFE.\nThis is not allowed.\n");
#endif

#if defined(GFM_RPROCESS) && !defined(GFM_CHEMTAGS)
  mpi_terminate("Code was compiled with GFM_RPROCESS, but not with GFM_CHEMTAGS.\nThis is not allowed.\n");
#endif

#if defined(GFM_DUST) && defined(GFM_SPLITFE)
  mpi_terminate("Code was compiled with GFM_SPLITFE, and GFM_DUST.\nThis is not supported yet.\n");
#endif

#if defined(GFM_CONST_IMF) && defined(GFM_VARIABLE_IMF)
  mpi_terminate("Code was compiled with GFM_CONST_IMF and GFM_VARIABLE_IMF together.\nThis is not allowed.\n");
#endif

#if defined(GFM_WINDS_VARIABLE) && !defined(FOF)
  mpi_terminate("Code was compiled with GFM_WINDS_VARIABLE, but not with FOF.\nThis is not allowed.\n");
#endif

#if(defined(GFM_BIPOLAR_WINDS) && (GFM_BIPOLAR_WINDS == 1)) && !defined(FOF)
  mpi_terminate("Code was compiled with GFM_BIPOLAR_WINDS == 1, but not with FOF.\nThis is not allowed.\n");
#endif

#if defined(GFM_WINDS_SAVE_PARTTYPE) && !defined(GFM_WINDS)
  mpi_terminate("GFM_WINDS_SAVE_PARTTYPE requires GFM_WINDS.\n");
#endif

#if defined(GFM_DISCRETE_ENRICHMENT) && (!defined(GFM) || !defined(GFM_STELLAR_EVOLUTION))
  mpi_terminate("GFM_DISCRETE_ENRICHMENT requires both GFM and GFM_STELLAR_EVOLUTION.\n");
#endif

#if defined(BH_BUBBLES) && !defined(BLACK_HOLES)
  mpi_terminate("Code was compiled with BH_BUBBLES, but not with BLACK_HOLES.\nThis is not allowed.\n");
#endif

#if defined(REPOSITION_ON_POTMIN) && !defined(EVALPOTENTIAL)
  mpi_terminate("REPOSITION_ON_POTMIN requires EVALPOTENTIAL.\n");
#endif

#if defined(BH_BASED_CGM_ZOOM) && (!defined(BLACK_HOLES) || !defined(REPOSITION_ON_POTMIN) || !defined(BH_DO_NOT_PREVENT_MERGERS))
  mpi_terminate("BH_BASED_CGM_ZOOM: Expect also BLACK_HOLES, REPOSITION_ON_POTMIN, and BH_DO_NOT_PREVENT_MERGERS.\n");
#endif

#if defined(VORONOI_FREQUENT_IMAGES) && defined(IMAGES_FOREACHSNAPSHOT)
  mpi_terminate("VORONOI_FREQUENT_IMAGES and IMAGES_FOREACHSNAPSHOT can't be used together");
#endif

#if defined(METALS) && !defined(USE_SFR) && !defined(ADDBACKGROUNDGRID)
  mpi_terminate("Code was compiled with METALS, but not with SFR.\nThis is not allowed.\n");
#endif

#if defined(MIN_METALLICITY_ON_STARTUP) && !defined(METALS)
  mpi_terminate("Code was compiled with MIN_METALLICITY_ON_STARTUP but not METALS.\nThis is not allowed.\n");
#endif

#if defined(TIMEDEPDE) && !defined(DARKENERGY)
  mpi_terminate("Code was compiled with TIMEDEPDE, but not with DARKENERGY.\nThis is not allowed.\n");
#endif

#if defined(OTVET) && defined(GFM_COOLING_METAL)
  mpi_terminate("Code was compiled with OTVET and GFM_COOLING_METAL together.\n This is not supported yet.\n");
#endif

#ifdef OTVET_FIXTIMESTEP
  if(All.MaxSizeTimestep != All.MinSizeTimestep)
    {
      mpi_terminate("OTVET_FIXTIMESTEP requires All.MaxSizeTimestep= All.MinSizeTimestep. Fix it in parameterfile. Stopping...");
    }
#endif

#if defined(OTVET_SCATTER_SOURCE) && defined(EDDINGTON_TENSOR_STARS) && defined(GFM)
  mpi_terminate("Code was compiled with OTVET_SCATTER_SOURCE for stars and GFM together.\n This is not supported yet.\n");
#endif

//#if defined(FM_RADIATION_FEEDBACK) && (!defined(FM_STAR_FEEDBACK) || !defined(GFM_STELLAR_EVOLUTION))
#if defined(FM_RADIATION_FEEDBACK) && (!defined(GFM_STELLAR_EVOLUTION))
  mpi_terminate(
      "Code was compiled with FM_RADIATION_FEEDBACK that requires FM_STAR_FEEDBACK and GFM_STELLAR_EVOLUTION active.\n Please "
      "activate both. Stopping.\n");
#endif

#ifdef MODIFIED_EOS
  check_modified_eos_parameters();
#endif

#if defined(BLACK_HOLES) && (!defined(BH_RELATIVE_NGB_DEVIATION))
  if(All.DesNumNgbBlackHole == All.MaxNumNgbDeviation)
    {
      mpi_terminate(
          "Code was compiled without BH_RELATIVE_NGB_DEVIATIONS and hence requires All.DesNumNgbBlackHole > All.MaxNumNgbDeviation.\n "
          "Fix choices in parameterfile or the code will crash. Stopping.\n");
    }
#endif

#undef REAL
#undef STRING
#undef INT
}

/*! \brief This function checks the consistency of the input parameters.
 *
 *  If you encounter some possible misuse and a corresponding error message
 *  that is hard to interpret, a check should be placed in this function with
 *  a terminate statement and a clear explanation why this does not work.
 *
 *  \return void
 */
void check_parameters()
{
  /* check whether time max is larger than max timestep */
  if(All.TimeMax - All.TimeBegin <= All.MaxSizeTimestep)
    {
      printf("PARAMETERS: check_parameters: TimeBegin = %g, TimeMax = %g, MaxSizeTimestep = %g \n", All.TimeBegin, All.TimeMax,
             All.MaxSizeTimestep);
      terminate(
          "check_parameters: Your total runtime is smaller or equal than the maximum allowed timestep! Choose an appropriate value "
          "for MaxSizeTimestep < TimeMax-TimeBegin! \n");
    }
}

/*! \brief This function reads a table with a list of desired output times.
 *
 *  The table does not have to be ordered in any way, but may not contain more
 *  than MAXLEN_OUTPUTLIST entries.
 *
 *  \param[in] fname The file name of the outputlist.
 *
 *  \return 0: success  1: unable to open file.
 */
int read_outputlist(char *fname)
{
  FILE *fd;
  int count, flag;
  char buf[512];

  if(!(fd = fopen(fname, "r")))
    {
      printf("can't read output list in file '%s'\n", fname);
      return 1;
    }

  All.OutputListLength = 0;

  while(1)
    {
      if(fgets(buf, 500, fd) != buf)
        break;

      count = sscanf(buf, " %lg %d ", &All.OutputListTimes[All.OutputListLength], &flag);

      if(count == 1)
        flag = 1;

      if(count == 1 || count == 2)
        {
          if(All.OutputListLength >= MAXLEN_OUTPUTLIST)
            terminate("\ntoo many entries in output list. You should increase MAXLEN_OUTPUTLIST=%d.\n", (int)MAXLEN_OUTPUTLIST);

          All.OutputListFlag[All.OutputListLength] = flag;
          All.OutputListLength++;
        }
    }

  fclose(fd);

  printf("\nBEGRUN: found %d times in output-list.\n", All.OutputListLength);

  return 0;
}
