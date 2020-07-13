#!/bin/bash
## shell script for code testing

#################################
## perform predefined examples ##
#################################

## Number of cores to compile and run the problem on
## use NUMBER_OF_TASKS=1 for 1d test problems!
NUMBER_OF_TASKS=1
NUMBER_OF_COMPILERS=8
PLOT=False # create plots for runs? True/False
PYTHON=python3

## choose your tests
TESTS=""
## available 1d test cases
TESTS+="wave_1d "
TESTS+="shocktube/shocktube_1d "
TESTS+="shocktube/shocktube_sod_1d "
TESTS+="interacting_blastwaves_1d "
TESTS+="polytrope_spherical_1d "
TESTS+="mhd_shocktube_1d "

## available 2d test cases
TESTS+="Gresho_2d "
TESTS+="Noh_2d/Noh_2d "
TESTS+="Noh_2d/Noh_refinement_2d "
TESTS+="Yee_2d "
TESTS+="khi_hydro_2d "
TESTS+="khi_mhd_2d "
TESTS+="current_sheet_2d "

## 3d test cases
TESTS+="Noh_3d/Noh_3d "
TESTS+="Noh_3d/Noh_refinement_3d "
#TESTS+="cosmo_box/gravity_only_3d "
#TESTS+="cosmo_box/star_formation_3d "
#TESTS+="cosmo_zoom_gravity_only_3d "
#TESTS+="isolated_galaxy_collisionless_3d "
#TESTS+="galaxy_merger_star_formation_3d "

## AMR test cases
#TESTS+="AMR/shocktube_2d "

## loop over all tests
for TEST in $TESTS
do
  DIR="examples/$TEST"
  RUNDIR="run/examples/$TEST"

  ## clean up
  rm -rf run/

  ## create run directory
  mkdir -p "$RUNDIR/"

  ## copy setup to run directory
  cp -r "$DIR"/* "$RUNDIR/"

  ## create ICs in run directory
  echo "$DIR/"
  $PYTHON "$DIR/create.py" "$RUNDIR/"
  ((return_value=$?))    ## get return value
  if [ $return_value != 0 ]; then    ## check return value
    printf 'ERROR: test.sh:\t%s\t python create.py failed!\n' "$DIR/"
    exit $return_value
  fi

  ## compile Arepo
  make -j "$NUMBER_OF_COMPILERS" CONFIG="$RUNDIR/Config.sh" BUILD_DIR="$RUNDIR/build" EXEC="$RUNDIR/Arepo"
  ((return_value=$?))    ## get return value
  if [ $return_value != 0 ]; then    ## check return value
    printf 'ERROR: test.sh:\t%s\t make failed!\n' "$DIR/"
    exit $return_value
  fi

  ## change to RUNDIR in subshell and execute test simulation
  (cd "$RUNDIR/" && mpiexec -n "$NUMBER_OF_TASKS" ./Arepo ./param.txt)
  ((return_value=$?))    ## get return value
  if [ $return_value != 0 ]; then    ## check return value
    printf 'ERROR: test.sh:\t%s\t execution failed!\n' "$DIR/"
    exit $return_value
  fi

  ## check result in example directory, this also creates some check plots
  $PYTHON "$DIR/check.py" "$RUNDIR/" "$PLOT"
  ((return_value=$?))    ## get return value
  if [ $return_value != 0 ]; then    ## check return value
    echo 'ERROR: test.sh: test failed!'
    exit $return_value
  else
    printf 'test.sh:\t%s\t test passed!\n' "$DIR/"
    echo 'cleaning up...'
  fi

  ## clean up
  rm -rf ./run
done

echo
echo
echo 'Tests'
echo $TESTS
echo 'passed!'

exit $return_value
