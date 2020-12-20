scriptsdir=../../Analysis/

tag=system
mintimestep=900000
rc=6.6667

for system in $(ls -1 | grep system_)
do
    mintimestep=900000
    if [[ ! -d ${system} ]];
    then
        continue
    fi
    echo $system
    check=$(echo $system | grep IHI | wc -l)
    if [[ $check -eq 1 ]];
    then
	mintimestep=1900000
	echo ${mintimestep}
    fi	

    for resolution in 30kb ;
    do 

        replicas="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49"
        check=$(ls -1 ./${system}/contact_map/distances_model_a_thaliana_replica_*_below_6_6667sigma_at_3kb_txt.tar.gz 2> /dev/null | wc -l | awk '{print $1}')
        echo $check
        if [[ ${check} -ge 49 ]];
        then
            replicas="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50"
        fi
        echo $replicas

	for replica in ${replicas} ;
	do
	    check=$(tail -1 $system/replica_${replica}/out* 2> /dev/null | grep Total | wc -l)
	    if [[ $check -eq 0 ]];
	    then
		echo "Simulation in replica $replica didn't finish!"
		continue
	    fi
            if [[ -e ./${system}/contact_map/norm_contacts_vs_gendist_a_thaliana_at_6.6667sigma_at_${resolution}_replica_${replica}_${system}.txt ]];		
            then
		ls -lrth ./${system}/contact_map/norm_contacts_vs_gendist_a_thaliana_at_6.6667sigma_at_${resolution}_replica_${replica}_${system}.txt
		rm -i jobscript_contact_${system}_${resolution}_replica_${replica}.cmd 2> /dev/null
		echo "Analysis in replica $replica is DONE!"
		continue
            fi
	    if [[ -e jobscript_contact_${system}_${resolution}_replica_${replica}.cmd ]];
            then
                echo "Analysis in replica $replica is ongoing!"
                continue
            fi

	    echo "Doing calculation in $replica"

	    script=00_contact_maps_genomewide.sh
	    
	    sed -e "s/XXXmintimestepXXX/${mintimestep}/g" -e "s/XXXresolutionXXX/${resolution}/g" -e "s/XXXsystemXXX/${system}/g" -e "s/XXXrcXXX/${rc}/g" -e "s/XXXreplicaXXX/${replica}/g" YYYscriptsdirYYY/${script} | bash

	done
    done
done # Close $system
