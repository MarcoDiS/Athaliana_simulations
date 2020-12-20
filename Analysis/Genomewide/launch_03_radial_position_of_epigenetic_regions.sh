scriptsdir=../../Analysis/
script=03_radial_position_of_epigenetic_regions.sh

dirs=$(ls -1 | grep system_N | grep -v log)

for dir in ${dirs} ;
do

    if [[ ! -d ${dir} ]];
    then
	continue
    fi

    if [[ ! -d ${dir}/contact_map/input_files ]];
    then
	continue
    fi

    cd ${dir}/contact_map/input_files

    #Telomere 1024
    #NORs 5058
    #Undetermined 10212
    #Polycomb_like 11940
    #Heterochromatin 12110
    #Active 44158
    for state in All NORs Active Telomere NORs Undetermined Polycomb_like Heterochromatin ;
    do
	echo $dir $state

	if [[ $state == "All" ]];
	then
            # Prepare DATA_FILES_INPUT.DAT from 1MLN up to 20MLN
	    ls -1 *_???????.XYZ  | grep XYZ >  DATA_FILES_INPUT.DAT    
	    ls -1 *_1???????.XYZ | grep XYZ >> DATA_FILES_INPUT.DAT    
	    ls -1 *_20000000.XYZ | grep XYZ >> DATA_FILES_INPUT.DAT    
	fi

	echo "Dir "$PWD
	echo "State "$state
	echo "Script "${script}

	bash ${scriptsdir}/${script} ${state}

    done

    cd ../../../

done
