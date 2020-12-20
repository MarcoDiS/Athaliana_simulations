for replica in $(seq 1 1 50);
do
    for radius in 2.2 ;
    do
	if [[ -d replica_${replica}_radius_${radius} ]];
	then
	    continue
	fi	
	echo replica_${replica}_radius_${radius}
    
	sed -e "s/XXXradiusXXX/${radius}/g" -e "s/XXXreplicaXXX/${replica}/g" jobscript_4_inflate_NADs_particles.cmd > _tmp_${replica}_radius_${radius}.cmd

	mnsubmit _tmp_${replica}_radius_${radius}.cmd | awk -v r=${replica} '{print "replica_"r,$NF}' >> running_jobs.log
    done
done
