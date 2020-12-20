for replica in $(seq 1 1 50);
do
    if [[ ! -e ../3_nucleolus_formation_vshape_from_yaxis_pulling_case/replica_${replica}/3_nucleolus_formation_part2_replica_${replica}.txt ]];
    then
	continue
    fi

    for radius in 2.2 ; # 0.64 of the nuclear volume expected for Random CLose Packing (RCP)
    do
	if [[ -d replica_${replica}_radius_${radius} ]];
	then
	    continue
	fi	
	echo replica_${replica}_radius_${radius}    
	
	bash 4_inflate_NADs_particles_vshape_case.sh ${replica} ${radius}
    done
done
