scriptsdir=../../Analysis/

dirs="system_NAD_ATTR1.0andREPU1.0_AC_ATTR0.20andNEUT_FH_ATTR0.50andNEUT_CH_NEUTandREPU0.000390625 system_NAD_ATTR1.0andREPU1.0_AC_ATTR0.20andNEUT_FH_NEUTandNEUT_CH_NEUTandREPU0.000390625"

for wdir in $dirs ;
do

    echo ${wdir}

    for cluster in $(grep Plot ${scriptsdir}/Utilities/H3K27me3_clusters.txt | awk '{print $3}' | grep 164);
    do
	strcluster=$(echo $cluster | sed "s/:/_/g")

	sed -e "s/XXXwdirXXX/${wdir}/g" -e "s,XXXclusterXXX,${cluster},g" ../scripts_and_programs/08_contact_maps_H3K27me3_geneclusters.sh | bash

    done
done
