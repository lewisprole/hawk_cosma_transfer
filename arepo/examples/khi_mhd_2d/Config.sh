# examples/khi_mhd_2d
NTYPES=6
NSOFTTYPES=2
PERIODIC
TWODIMS
GAMMA=1.66666666666666
VORONOI
REGULARIZE_MESH_CM_DRIFT
REGULARIZE_MESH_CM_DRIFT_USE_SOUNDSPEED
REGULARIZE_MESH_FACE_ANGLE
#TREE_BASED_TIMESTEPS # off
DOUBLEPRECISION=1
OUTPUT_IN_DOUBLEPRECISION # snapshot files will be written in double precision
INPUT_IN_DOUBLEPRECISION  # initial conditions are in double precision

# generic
VORONOI_DYNAMIC_UPDATE
HAVE_HDF5
DEBUG
OUTPUT_VORTICITY

# PASSIVE_SCALARS=1

LONG_X=0.5792653228386294

FORCE_EQUAL_TIMESTEPS

# LONG_Y=2.0

VORONOI_STATIC_MESH

MESHRELAX_DENSITY_IN_INPUT

OUTPUT_CENTER_OF_MASS

MHD
RIEMANN_HLLD

TETRA_INDEX_IN_FACE 
# MHD_POWELL
