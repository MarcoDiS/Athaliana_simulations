dir=${PWD}

for wdir in $(ls -1 | grep system_NAD_ATTR1.0andREPU1.0_ | grep -v "sh\|log\|cmd")
do
    echo $wdir
    if [[ ! -d ${wdir}/contact_map/ ]];
    then
        continue
    fi
    if [[ ! -e ${wdir}/contact_map/contact_matrix_a_thaliana_at_6.6667sigma_at_30kb_${wdir}.tab ]];
    then                                                                                                          
	echo "Contact matrix not present!"
	continue
    fi 

    script=02_compartment_strength_analysis_on_models_singlechromosome.sh

    sed -e "s,XXXwdirXXX,${wdir},g" -e "s,XXXdirXXX,${dir},g" ../../Analysis/SingleChromosome/${script} | bash

done # Close cycle over $dir
