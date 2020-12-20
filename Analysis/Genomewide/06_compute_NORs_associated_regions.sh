scriptsdir=../../Analysis/
input=${scriptsdir}/contact_maps/

wdir=contact_map
format=tab
rc=6.6667

nreplicas=50

for dir in XXXwdirXXX ;
do
    echo $dir
    echo $wdir
    cd $dir  2> /dev/null
    cd $wdir 2> /dev/null
    alltmpfile=_tmp_all_distances
    NORstmpfile=_tmp_NORs_bins
    noNORstmpfile=_tmp_noNORs_bins
    outfile=contacts_with_NORs_at_${rc}sigma_at_3kb_${dir}.tab

    gzfile=distances_model_a_thaliana_replica_${replica}_below_6_6667sigma_at_3kb_txt.tar.gz
    if [[ ! -e ${alltmpfile} ]];
    then

	zcat distances_model_a_thaliana_replica_*_below_6_6667sigma_at_3kb_txt.tar.gz > ${alltmpfile}
    fi
    ls -lrth ${alltmpfile}    
    
    grep NORs ${scriptsdir}/Utilities/chromatin_states_haploid_genome_with_locations.txt | awk '{print "NORs",$2-1}' > ${NORstmpfile}
    grep -v NORs ${scriptsdir}/Utilities/chromatin_states_haploid_genome_with_locations.txt | awk '{print "noNORs",$1,$2-1}' > ${noNORstmpfile}
    
    awk '{h[$1]+=$2}END{for(i in h) print i,h[i]}' <(awk '{h[$1]=0}END{for(i in h) print i,h[i]}' <(grep chr ${scriptsdir}/Utilities/chromatin_states_haploid_genome_with_locations.txt)) <( awk -v rc=$rc '{if($1=="NORs"){NORs[$2]=$1; next}; if($1=="noNORs"){noNORs[$3]=$2; next}; {if(($2 in NORs) && ($3 in NORs)){next}; if($4>rc*rc){next}; if($2 in NORs){print noNORs[$3],1}; if($3 in NORs){print noNORs[$2],1}}}' ${NORstmpfile} ${noNORstmpfile} ${alltmpfile}) > ${outfile}
    rm ${NORstmpfile} ${noNORstmpfile} 
    cd .. # Exit ${wdir}
    cd .. # Exit ${dir}
    
done
echo "Analysis DONE!"
