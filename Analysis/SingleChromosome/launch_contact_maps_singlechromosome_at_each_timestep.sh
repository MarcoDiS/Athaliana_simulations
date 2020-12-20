scriptsdir=../../Analysis/

rc=6.6667

for system in $(ls -1 | system_NAD_ATTR1.0andREPU1.0 | grep -v 'log|sh|txt') ;
do
    for resolution in 30kb ;
    do 

	echo ${system} ${resolution}

	script=contact_maps_singlechromosome_at_each_timestep.sh

	sed -e "s,XXXresolutionXXX,${resolution},g" -e "s,XXXwdirXXX,${wdir},g" -e "s,XXXrcXXX,${rc},g" ${scriptsdir}/${script} | bash
	
    done
done

