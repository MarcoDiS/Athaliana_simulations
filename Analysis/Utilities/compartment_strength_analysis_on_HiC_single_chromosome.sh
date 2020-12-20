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

for chr in chr4 ; #chr1 chr2 chr3 chr4 chr5 ;
do
    (
	c=$(echo $chr | sed "s/chr//g")
	statefile=_HiC_states_chr${c}
	matrixfile=_HiC_matrix_chr${c}

	#inmatrix=/home/devel/mstefano/cumulative_normalized_HiC_matrices_A_thaliana/cumulative_norm_HiC_matrices_${c}_dna_ObsOverExp_at_30kb.tab
	inmatrix=/home/devel/mstefano/cumulative_normalized_HiC_matrices_A_thaliana/cumulative_norm_HiC_matrices_${c}_dna_ObsOverExp_at_30kb_using_median.tab
	
	name=$(echo ${inmatrix} | sed -e "s,/, ,g" -e "s/\.tab//g" | awk '{print $NF}')
	#outfile=compartment_strength_cumulative_norm_HiC_matrices_${c}_dna_at_30kb.txt
	outfile=compartment_strength_cumulative_norm_HiC_matrices_${c}_dna_at_30kb_using_median.txt
	echo $instates $inmatrix $outfile
	
	if [[ ! -s ${outfile} ]];
	then
	    echo "#Bin state CS(bin)" | awk '{printf("%-9s\t%-8s\t%s\n", $1,$2,$3)}' > ${outfile}
	fi
	awk '{print $1"_"$2,$3"_"$4,$5}' <(cat ${inmatrix} | grep -v "MASKED\|#") > ${matrixfile}

	sed "s/:/ /g" ${instates} | grep -v "#" | awk '{print $1"_"$3,$4}' > ${statefile}
	offset=$(grep ${chr} ${statefile} | sed "s/_/ /g" | head -1 | awk '{print $2}')
	echo "Offset "$offset

	for bin in $( awk '{print $1}' ${statefile} | grep ${chr} ); 
	do

	    state=$(awk -v b=${bin} '{if(b==$1) print $2}' ${statefile})
	    
	    bin1=$(echo $bin | sed "s/_/ /g" | awk -v o=${offset} '{print $1"_"$2-o}')

	    check=$(grep -w ${bin1} ${outfile} | wc -l | awk '{print $1}')
	    if [[ $check -gt 0 ]];
	    then
		continue
	    else
		echo $bin1
	    
		CS=$(awk '{if(NF==2){state[$1]=$2}else{print $1,state[$1],$2,state[$2],$3}}' <( awk -v o=${offset} '{print $1"_"$2-o,$3}' <( sed "s/_/ /g" ${statefile} )) <( cat ${matrixfile} | grep -w ${bin1}) | awk '{if($2==$4){num+=$5;cnum++}; den+=$5;cden++}END{print (num/cnum)/(den/cden),cnum,cden}' | awk '{print $1}')
		echo ${bin1} ${state} ${CS} | awk '{printf("%-9s\t%-8s\t%s\n", $1,$2,$3)}' >> ${outfile}
		tail -1 ${outfile}
	    fi
	done

    ) > output_${chr}_single_chromosome.log &
done