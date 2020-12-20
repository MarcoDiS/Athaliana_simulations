scriptsdir=../../Analysis/

mintimestep=900000
rc=6.6667 # 200nm
script=00_contact_maps_singlechromosome.sh


for system in $(ls -1 | grep system_NAD_ATTR1.0andREPU1.0_)
do
    if [[ ! -d ${system} ]];
    then
	continue
    fi
    echo $system

    for resolution in 30kb ;
    do 
	replicas="1 2 3 4 5 6 7 8 9"
	check=$(ls -1 ./${system}/contact_map/distances_model_a_thaliana_replica_*_below_6_6667sigma_at_3kb_txt.tar.gz 2> /dev/null | wc -l | awk '{print $1}') 
	echo $check
	if [[ ${check} -ge 9 ]];
	then
	    replicas="1 2 3 4 5 6 7 8 9 10"
	fi
	echo $replicas

	for replica in ${replicas} ;
	do

	    echo $system ${resolution} ${replica}
	    
	    sed -e "s/XXXmintimestepXXX/${mintimestep}/g" -e "s/XXXresolutionXXX/${resolution}/g" -e "s/XXXsystemXXX/${system}/g" -e "s/XXXrcXXX/${rc}/g" -e "s/XXXreplicaXXX/${replica}/g" ${scriptsdir}/${script} | bash
	    
	done # Close cycle over $replica
    done # Close cycle over $resolution
done # Close cycle over $system
