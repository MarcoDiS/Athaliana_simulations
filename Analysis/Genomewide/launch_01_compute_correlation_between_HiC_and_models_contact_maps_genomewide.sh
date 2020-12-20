scriptsdir=../../Analysis/

dirs=$(ls -1 | grep system_NAD)

for d in ${dirs} ;
do
    if [[ ! -d ${d} ]];
    then
        continue
    fi

    for replica in 0 ;
    do

	output=output_compute_correlation_between_contact_maps.log

	if [[ ! -e correlations_of_genomewide_contact_maps_in_genomewide_models_replica_${replica}_${d}.log ]];
	then
	    script=01_compute_correlation_between_HiC_and_models_contact_maps_genomewide_genomewide.sh
	    echo ${script}
	    bash ${scriptsdir}/${script} ${replica} ${d}
	fi
	
	if [[ ! -e correlations_of_intra_chromosome_contact_maps_in_genomewide_models_replica_${replica}_${d}.log ]];
	then
	    script=01_compute_correlation_between_HiC_and_models_contact_maps_genomewide_intra_chromosome.sh
	    echo ${script}
	    bash ${scriptsdir}/${script} ${replica} ${d}
	fi
	
	if [[ ! -e correlations_of_genomewide_contacts_vs_L_in_genome_wide_models_replica_${replica}_${d}.log ]];
	then
	    script=01_compute_correlation_between_HiC_and_models_contact_maps_genomewide_contacts_vs_L.sh
	    echo ${script}
	    bash ${scriptsdir}/${script} ${replica} ${d}
	fi
	echo "Analysis DONE!"

    done

done
