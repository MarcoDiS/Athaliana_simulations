
dirs="system_NAD_ATTR1.0andREPU1.0_AC_ATTR0.20andNEUT_FH_ATTR0.50andNEUT_CH_NEUTandREPU0.000390625"

for wdir in $dirs ;
do
    echo ${wdir}

    sed "s/XXXwdirXXX/${wdir}/g" ../../Analysis/06_compute_NORs_associated_regions.sh | bash

done
