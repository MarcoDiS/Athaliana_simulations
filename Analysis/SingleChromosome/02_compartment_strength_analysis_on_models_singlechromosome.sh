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

OutfileOverall=${PWD}/pvalues_parametrization_Wilcoxon_overall.txt

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

name=$(echo ${inmatrix} | sed -e "s,/, ,g" -e "s/\.tab//g" | awk '{print $NF"_single_chromosome"}')
echo $instates $inmatrix $outfile

if [[ ! -s ${outfile} ]];
then
    echo "#Bin state CS(bin)" | awk '{printf("%-9s\t%-8s\t%s\n", $1,$2,$3)}' > ${outfile}

    for chr in chr4 ;
    do
	c=$(echo $chr | sed "s/chr//g")
	statefile=_model_states_chr${c}
	matrixfile=_model_matrix_chr${c}
	
        # Annotate MASKED columns
	paste <(cat ${inmatrix} | grep -v "#") <( cat ../../Analysis/Utilities/*/cumulative_norm_HiC_matrices_${c}_dna_ObsOverExp_at_30kb.tab | grep -v "#" | awk '{print $NF}') | awk '{print $1,$2,$3,$4,$5,$7}' > _tmp_masked_matrix

        echo "#Gendist avg_contacts for genome-wide matrix"
        indecay=median_number_of_contacts_vs_gendist_models_${chr}_at_30kb.txt

        rm -fvr ${indecay}
        tmpfile=_tmp_sorted_${chr}
        awk -v chr=${chr} -v res=${resolution} '{if($1==chr && $3==chr){d=sqrt(($2-$4)*($2-$4)); print d*res,$5}else{print "Inter",$5}}' <( cat _tmp_masked_matrix | grep -v "MASKED" ) | sort -k 2n,2n > ${tmpfile}
        for d in $(awk '{print $1}' ${tmpfile} | uniq | sort -k 1n | uniq);
        do
            awk -v d=$d 'BEGIN{cnt=0}{if($1==d){cnt++;v[cnt]=$2}}END{if(cnt%2==1){print d,v[int(cnt*0.5)+1]}else{print d,(v[int(cnt*0.5)+1]+v[int(cnt*0.5)])*0.5}}' ${tmpfile} >> ${indecay}
            tail -1 ${indecay}
        done

        # Compute ObsOverExp matrix
	outmatrix=${inmatrix%.tab}_ObsOverExp_using_median.tab
	echo $outmatrix
	
	echo "#Gendist avg_contacts for ${chr}"
	grep "#" ${inmatrix} > ${outmatrix}
	head ${indecay} ${inmatrix}
	awk -v res=${resolution} -v chr=${chr} '{if(NF==2){e[$1]=$2}else{if($1==chr && $3==chr){d=int(sqrt(($2-$4)*($2-$4))*res); print $1,$2,$3,$4,$5/e[d],$6};}}' ${indecay} _tmp_masked_matrix >> ${outmatrix}
	
	inmatrix=${inmatrix%.tab}_ObsOverExp_using_median.tab
	
        # Chromatin state per bin
	sed "s/:/ /g" ${instates} | grep -v "#" | awk '{print $1"_"$3,$4}' > ${statefile}
	offset=$(grep ${chr} ${statefile} | sed "s/_/ /g" | head -1 | awk '{print $2}')
	echo "Offset "$offset
	
	for bin in $(awk '{print $1}' ${statefile} | grep ${chr}); 
	do
	    check=$(grep -w ${bin} ${outfile} | wc -l | awk '{print $1}')
	    if [[ $check -gt 0 ]];
	    then
		echo "Bin ${bin} DONE!"
		continue
	    fi
	    echo $bin
	    state=$(sed "s/:/ /g" ${instates} | awk -v b=${bin} '{if(b==$1"_"$3)print $4}')
	    
	    bin1=$(echo $bin | sed "s/_/ /g" | awk -v o=${offset} '{print $1"_"$2-o}')
	    echo $bin1
	    CS=$(awk '{if(NF==2){state[$1]=$2}else{print $1,state[$1],$2,state[$2],$3}}' <( awk -v o=${offset} '{print $1"_"$2-o,$3}' <( sed "s/_/ /g" ${statefile} )) <( awk '{if($1==$3) print $1"_"$2,$3"_"$4,$5}' <(cat ${inmatrix} | grep -v MASKED) | grep -w ${bin1}) | awk '{if($2==$4){num+=$5;cnum++}; den+=$5;cden++}END{print (num/cnum)/(den/cden),cnum,cden}' | awk '{print $1}')
	    
	    echo ${bin} ${state} ${CS} | awk '{printf("%-9s\t%-8s\t%s\n", $1,$2,$3)}' >> ${outfile}
	    tail -1 ${outfile}
	    
	done
    done
fi

rm _*


