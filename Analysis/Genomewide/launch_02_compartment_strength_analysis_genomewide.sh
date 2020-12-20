scriptsdir=../../Analysis/
dir=${PWD}

# Genome-wide
wdirs=$(ls -1 | grep system_N | grep -v log)

for wdir in ${wdirs} ;
do
    if [[ ! -d ${wdir} ]];
    then
	echo "${wdir} doesn't exist"
	continue
    fi
    if [[ ! -d ${wdir}/contact_map/ ]];
    then
	echo "${wdir}/contact_map doesn't exist"
	continue
    fi
    if [[ -e ${wdir}/contact_map/compartment_strength_contact_matrix_a_thaliana_at_6.6667sigma_at_30kb_${wdir}.txt ]];
    then
	echo "Analysis DONE!"
	continue
    fi
    echo $wdir
    script=compartment_strength_analysis_genomewide.sh

    sed -e "s,XXXwdirXXX,${wdir},g" -e "s,XXXdirXXX,${dir},g" ${scriptsdir}/${script} | bash

done # Close cycle over $dir
