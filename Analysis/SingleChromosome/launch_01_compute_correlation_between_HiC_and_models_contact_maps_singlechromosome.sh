scriptsdir=../../Analysis/

rc=6.6667

for system in $(ls -1 | grep system_NAD_ATTR1.0andREPU1.0_) ;
do
    for resolution in 30kb ;
    do 

	echo ${system} ${resolution}

	script=01_compute_correlation_between_HiC_and_models_contact_maps_singlechromosome.sh

	sed -e "s/XXXresolutionXXX/${resolution}/g" -e "s/XXXsystemXXX/${system}/g" -e "s/XXXrcXXX/${rc}/g" ${scriptsdir}/${script} | bash
	
    done # Close cycle over $resolution
done # Close cycle over $system

