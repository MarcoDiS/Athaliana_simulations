
dirs="system_NAD_ATTR1.0andREPU1.0_AC_ATTR0.20andNEUT_FH_ATTR0.50andNEUT_CH_NEUTandREPU0.000390625"

for wdir in $dirs ;
do
    for tag in centromeres ;
    do

	echo ${wdir} ${tag}

	script=07_position_of_interesting_genomic_regions.sh
	sed "s/XXXwdirXXX/${wdir}/g" ../../Analysis/Genomewide/${script} | bash
	
	script=07_${tag}_mutual_distances.sh	
	sed "s/XXXwdirXXX/${wdir}/g" ../../Analysis/Genomewide/${script} | bash

    done
done
