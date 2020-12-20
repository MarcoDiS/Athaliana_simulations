#!/bin/bash

# @ job_name         = rep_XXXreplicaXXX
# @ initialdir       = /scratch/devel/mstefano/2018_05_04_Project_A_thaliana_physics_models/WT_simulations/Initial_conformations/generate_initial_conformation_v_chromosomes/4_inflate_NADs_particles_linear_case
# @ output           = output.log
# @ error            = error.log
# @ total_tasks      = 1
# @ cpus_per_task    = 8
# @ wall_clock_limit = 6:00:00
# @ partition        = genB,main
# @ class            = 4dgenome
# @ requeue          = 1

module purge
module load openmpi


bash 4_inflate_NADs_particles_linear_case.sh XXXreplicaXXX XXXradiusXXX
