#!/bin/bash

#SBATCH --job-name=collapsE                     # Job name
#SBATCH --mail-type=END,FAIL                     # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=prolel@cardiff.ac.uk    # Where to send mail
#SBATCH -p compute				#compute is the main queue, dev is the developement queue, dev requires ntasks=40, time=01:00:00
#SBATCH --ntasks=80                              # Number of MPI ranks
#SBATCH --ntasks-per-node=40                     # How many tasks on each node
#SBATCH --distribution=cyclic:cyclic             # Distribute tasks cyclically on nodes and sockets
#SBATCH -e error.txt
#SBATCH -o output.txt
#SBATCH --mem-per-cpu=4000mb                      # Memory per processor
#SBATCH --time=72:00:00                          # Time limit hrs:min:sec
#SBATCH --output=collapse%j.log                 # Standard output and error log
#SBATCH --exclusive
#SBATCH --account=scw1211

# change to current working directory
cd $SLURM_SUBMIT_DIR



# set memory limits
ulimit -l unlimited

# load any modules that moght be needed

module load intel
module load mpi/intel.4.1.0
module load gmp/5.0.2
module load gsl/1.15

export I_MPI_FABRICS=shm:ofa

# run the job
mpirun ./Arepo_POPIII_MHDvariable_sink_16 inArepo_popIII.param 0 > arepo_out.txt

