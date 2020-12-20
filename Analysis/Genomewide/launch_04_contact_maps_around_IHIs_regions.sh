scriptsdir=../../Analysis/

dirs="system_NAD_ATTR1.0andREPU1.0_AC_ATTR0.20andNEUT_FH_ATTR0.50andNEUT_CH_NEUTandREPU0.000390625_IHIs_springs_from_optimal_conformation system_NAD_ATTR1.0andREPU1.0_AC_ATTR0.20andNEUT_FH_ATTR0.50andNEUT_CH_NEUTandREPU0.000390625"


for wdir in $dirs ;
do

    echo ${wdir}

    for nIHI1 in $(seq 1 1 10);
    do
	for nIHI2 in $(seq $((${nIHI1}+1)) 1 10);
	do

	    script=04_contact_maps_around_IHIs_regions.sh
	    sed -e "s/XXXwdirXXX/${wdir}/g" -e "s/XXXnIHI1XXX/${nIHI1}/g" -e "s/XXXnIHI2XXX/${nIHI2}/g" ${scriptsdir}/${script} | bash
	   	    
	done
    done
done
