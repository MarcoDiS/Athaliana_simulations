scriptsdir=/scratch/devel/mstefano/2018_05_04_Project_A_thaliana_physics_models/WT_simulations/Initial_conformations/generate_initial_conformation_v_chromosomes/scripts_and_programs/

mintimestep=900000
rc=6.6667

templatefile=${scriptsdir}/jobscript_contact_maps_a_thaliana_genomewide_at_fixed_timestep.cmd

for wdir in system_NAD_ATTR1.0andREPU1.0_AC_ATTR0.20andNEUT_FH_ATTR0.50andNEUT_CH_NEUTandREPU0.000390625 ;
do
    echo $wdir

    for resolution in 30kb ;
    do 

	echo ${resolution}

	jobfile=jobscript_contact_${wdir}_${resolution}_at_fixed_timestep.cmd

	sed -e "s,YYYdirYYY,${PWD},g" -e "s,YYYscriptsdirYYY,${scriptsdir},g" -e "s/YYYresolutionYYY/${resolution}/g" -e "s/YYYwdirYYY/${wdir}/g" -e "s/YYYrcYYY/${rc}/g" ${templatefile} > ${jobfile}
	mnsubmit ${jobfile}
	
    done
done

