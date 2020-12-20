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
# /home/devel/mstefano/cumulative_normalized_HiC_matrices_A_thaliana/cumulative_norm_HiC_matrices_full_at_30kb.tab

instates=/scratch/devel/mstefano/2018_05_04_Project_A_thaliana_physics_models/WT_simulations/Initial_conformations/generate_initial_conformation_v_chromosomes/scripts_and_programs/chromatin_states_haploid_genome_with_locations_at_30kb.txt
#inmatrix=/home/devel/mstefano/cumulative_normalized_HiC_matrices_A_thaliana/cumulative_norm_HiC_matrices_full_ObsOverExp_at_30kb.tab
inmatrix=/home/devel/mstefano/cumulative_normalized_HiC_matrices_A_thaliana/cumulative_norm_HiC_matrices_full_ObsOverExp_at_30kb_using_median.tab
name=$(echo ${inmatrix} | sed -e "s,/, ,g" -e "s/\.tab//g" | awk '{print $NF}')
#outfile1=compartment_strength_cumulative_norm_HiC_matrices_full_at_30kb_only_intrachromosome.txt
#outfile2=compartment_strength_cumulative_norm_HiC_matrices_full_at_30kb_genomewide.txt
outfile1=compartment_strength_cumulative_norm_HiC_matrices_full_at_30kb_using_median_only_intrachromosome.txt
outfile2=compartment_strength_cumulative_norm_HiC_matrices_full_at_30kb_using_median_genomewide.txt
echo $instates $inmatrix $outfile

if [[ ! -s ${outfile1} ]];
then
    echo "#Bin state CS(bin)" | awk '{printf("%-9s\t%-8s\t%s\n", $1,$2,$3)}' > ${outfile1}
fi
if [[ ! -s ${outfile2} ]];
then
    echo "#Bin state CS(bin)" | awk '{printf("%-9s\t%-8s\t%s\n", $1,$2,$3)}' > ${outfile2}
fi

awk '{if($1==$3) print $1"_"$2,$3"_"$4,$5}' <(cat ${inmatrix} | grep -v "MASKED\|#") > _HiC_only_intrachromosome
awk '{print $1"_"$2,$3"_"$4,$5}' <(cat ${inmatrix} | grep -v "MASKED\|#") > _HiC_genomewide
sed "s/:/ /g" ${instates} | grep -v "#" | awk '{print $1"_"$3,$4}' > _HiC_states


for chr in chr1 chr2 chr3 chr4 chr5 ;
do
    (
	for bin in $( awk '{print $1}' _HiC_states | grep ${chr} ); 
	do
	    state=$(awk -v b=${bin} '{if(b==$1) print $2}' _HiC_states)

	    check=$(grep -w ${bin} ${outfile1} | wc -l | awk '{print $1}')
	    if [[ $check -gt 0 ]];
	    then
		continue
	    else
		echo $bin
	    
		CS=$(awk '{if(NF==2){state[$1]=$2}else{print $1,state[$1],$2,state[$2],$3}}' _HiC_states <( cat _HiC_only_intrachromosome | grep -w ${bin}) | awk '{if($2==$4){num+=$5;cnum++}; den+=$5;cden++}END{print (num/cnum)/(den/cden),cnum,cden}' | awk '{print $1}')
		echo ${bin} ${state} ${CS} | awk '{printf("%-9s\t%-8s\t%s\n", $1,$2,$3)}' >> ${outfile1}
		tail -1 ${outfile1}
	    fi

	    check=$(grep -w ${bin} ${outfile2} | wc -l | awk '{print $1}')
	    if [[ $check -gt 0 ]];
	    then
		continue
	    else
		echo $bin
		CS=$(awk '{if(NF==2){state[$1]=$2}else{print $1,state[$1],$2,state[$2],$3}}' _HiC_states <( cat _HiC_genomewide      | grep -w ${bin}  ) | awk '{if($2==$4){num+=$5;cnum++}; den+=$5;cden++}END{print (num/cnum)/(den/cden),cnum,cden}' | awk '{print $1}')
		echo ${bin} ${state} ${CS} | awk '{printf("%-9s\t%-8s\t%s\n", $1,$2,$3)}' >> ${outfile2}
		tail -1 ${outfile2}
	    fi

	done

    ) > output_${chr}.log &
done
