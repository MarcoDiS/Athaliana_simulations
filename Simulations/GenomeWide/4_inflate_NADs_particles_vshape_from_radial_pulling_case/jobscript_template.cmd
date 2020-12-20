#!/bin/bash

# @ job_name         = rep_XXXreplicaXXX
# @ initialdir       = XXXdirXXX
# @ output           = output_replica_XXXreplicaXXX.log
# @ error            = error.log
# @ total_tasks      = 1
# @ cpus_per_task    = 8
# @ wall_clock_limit = XXXtimeXXX:00:00
# @ partition        = genB,main
# @ class            = XXXclassXXX
# @ requeue          = 1

module purge
module load openmpi

time ( mpirun -np 8 ~/LAMMPS/lammps-31Mar17_parallel_version/src/lmp_mpi -in _tmp.lmp -log none -echo screen )
