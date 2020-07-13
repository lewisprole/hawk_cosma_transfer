# examples/AMR/KH_2d

NTYPES=6
PERIODIC
TWODIMS

AMR
#VORONOI

AMR_GRADIENTS
#AMR_CONNECTIONS

TVD_SLOPE_LIMITER
TVD_SLOPE_LIMITER_SUPERBEE

FORCE_EQUAL_TIMESTEPS

REFINEMENT_SPLIT_CELLS
REFINEMENT_MERGE_CELLS

NO_GAS_SELFGRAVITY

DOUBLEPRECISION=1
DOUBLEPRECISION_FFTW
OUTPUT_IN_DOUBLEPRECISION
INPUT_IN_DOUBLEPRECISION

OUTPUT_PRESSURE_GRADIENT
OUTPUT_DENSITY_GRADIENT
OUTPUT_VELOCITY_GRADIENT
OUTPUT_SURFACE_AREA
OUTPUT_PRESSURE

# generic
CHUNKING
HAVE_HDF5
DEBUG
