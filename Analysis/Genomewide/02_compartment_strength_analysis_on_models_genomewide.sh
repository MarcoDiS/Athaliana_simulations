# Given an assignment of bins to compartments, we define compartment
# strength first on a per-bin level. The compartment strength of bin i
# (CSi) is the average number of contacts it makes with other bins of
# the same compartment type in the observed-over-expected heat map,
# divided by the average number of contacts it makes with any bin in the
# observed-over-expected heat map. The compartment strength of the total
# data set is then <CSi>, where the average is taken over all bins,
# weighted equally. Note that this metric is independent of the
# orientation of the compartment profile, since the two compartments are
# treated symmetrically. If there is no compartmentalization, the metric
# is 1, whereas any pattern of compartmentalization yields a compartment
# strength greater than 1.
# State strength!!!

# chromatin_states_haploid_genome_with_locations_at_30kb.txt
dir=XXXdirXXX
wdir=XXXwdirXXX
pwd
cd ${dir} 
pwd
cd ${wdir}/contact_map/
pwd
resolution=30000

instates=../../Analysis/Utilities/chromatin_states_haploid_genome_with_locations_at_30kb.txt
inmatrix=$(ls contact_matrix_a_thaliana_at_6.6667sigma_at_30kb_s*tab | grep -v Obs 2> /dev/null)

name=$(echo ${inmatrix} | sed -e "s,/, ,g" -e "s/\.tab//g" | awk '{print $NF}')
outfile=compartment_strength_${name}.txt
echo $instates $inmatrix $outfile

if [[ ! -s ${outfile} ]];
then
    echo "#Bin state CS(bin)" | awk '{printf("%-9s\t%-8s\t%s\n", $1,$2,$3)}' > ${outfile}
fi

# Annotate MASKED columns
paste <(cat ${inmatrix} | grep -v "#") <( cat ../../Analysis/Utilities/*/cumulative_norm_HiC_matrices_full_ObsOverExp_at_30kb.tab | grep -v "#" | awk '{print $NF}') | awk '{print $1,$2,$3,$4,$5,$7}' > _tmp_masked_matrix

# Expected number of contacts
indecay=average_number_of_contacts_vs_gendist_models_${wdir}_at_30kb.txt
awk -v res=${resolution} '{if($1==$3){d=sqrt(($2-$4)*($2-$4))}else{d="Inter"}; h[d]+=$5; h2[d]+=$5*$5; cnt[d]++}END{for(i in h){avg=h[i]/cnt[i]; avg2=h2[i]/cnt[i]; stddev=sqrt(avg2-avg*avg); if(i=="Inter"){print i,avg,stddev}else{print (i)*res,avg,stddev}}}' <( cat _tmp_masked_matrix | grep -v "MASKED" ) | sort -k 1n > ${indecay}

# Compute ObsOverExp matrix
outmatrix=${inmatrix%.tab}_ObsOverExp.tab
echo $outmatrix

echo "#Gendist avg_contacts"
grep "#" ${inmatrix} > ${outmatrix}
head ${indecay} ${inmatrix}
awk -v res=${resolution} '{if(NF==3){e[$1]=$2}else{if($1==$3){d=int(sqrt(($2-$4)*($2-$4))*res); print $1,$2,$3,$4,$5/e[d],$6}else{d="Inter"; print $1,$2,$3,$4,$5/e[d],$6};}}' ${indecay} _tmp_masked_matrix >> ${outmatrix}

inmatrix=${inmatrix%.tab}_ObsOverExp.tab

awk '{if($1==$3) print $1"_"$2,$3"_"$4,$5}' <(cat ${inmatrix} | grep -v "MASKED\|#") > _model

# Chromatin state per bin
statefile=_model_states
sed "s/:/ /g" ${instates} | grep -v "#" | awk '{print $1"_"$3,$4}' > ${statefile}

for bin in $(awk '{print $1}' ${statefile}); 
do
    state=$(sed "s/:/ /g" ${instates} | awk -v b=${bin} '{if(b==$1"_"$3)print $4}')

    check=$(grep -w ${bin} ${outfile} | wc -l | awk '{print $1}')
    if [[ $check -gt 0 ]];
    then
	echo "Bin ${bin} DONE!"
        continue
    else
	echo $bin
	CS=$(awk '{if(NF==2){state[$1]=$2}else{print $1,state[$1],$2,state[$2],$3}}' ${statefile} <( cat _model | grep -w ${bin}) | awk '{if($2==$4){num+=$5;cnum++}; den+=$5;cden++}END{print (num/cnum)/(den/cden),cnum,cden}' | awk '{print $1}')
    
	echo ${bin} ${state} ${CS} | awk '{printf("%-9s\t%-8s\t%s\n", $1,$2,$3)}' >> ${outfile}
    fi
    tail -1 ${outfile}

done
rm _*
