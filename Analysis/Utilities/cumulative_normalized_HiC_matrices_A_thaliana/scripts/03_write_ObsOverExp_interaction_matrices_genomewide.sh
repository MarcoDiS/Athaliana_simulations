# Compute contacts vs L
for res in 30kb ; 
do
    echo $res
    resolution=$(echo ${res} | sed 's/kb/000/g')
    inmatrix=cumulative_norm_HiC_matrices_full_at_${res}.tab
    indecay=average_number_of_contacts_vs_gendist_a_thaliana_norm_full_${res}.txt

    # Genome-wide
    echo "#Gendist avg_contacts for genome-wide matrix"
    outmatrix=cumulative_norm_HiC_matrices_full_ObsOverExp_at_${res}.tab
    grep "#" ${inmatrix} > ${outmatrix}
    head ${indecay} ${inmatrix} 
    awk -v res=${resolution} '{if(NF==3){e[$1]=$2}else{if($1==$3){d=int(sqrt(($2-$4)*($2-$4))*res); print $1,$2,$3,$4,$5/e[d],$6}; if($1!=$3){d="Inter"; print $1,$2,$3,$4,$5/e[d],$6}}}' ${indecay} <( grep -v "#" ${inmatrix} ) >> ${outmatrix}

    # Per chromosome
    for chr in chr1 chr2 chr3 chr4 chr5 ;
    do
	c=$(echo $chr | sed "s/chr//g" | awk '{print $1"_dna"}')

	indecay=average_number_of_contacts_vs_gendist_a_thaliana_norm_${c}_at_${res}.txt

	outmatrix=cumulative_norm_HiC_matrices_${c}_ObsOverExp_at_${res}.tab
	
	echo "#Gendist avg_contacts for ${chr}"
	grep "#" ${inmatrix} > ${outmatrix}
	head ${indecay} ${inmatrix} 
	offset=$(awk -v chr=${chr} '{if($1==chr && $3==chr){print $0}}'  <( grep -v "#" ${inmatrix} ) | head -1 | awk '{print $2}')
	awk -v o=${offset} -v res=${resolution} -v chr=${chr} '{if(NF==3){e[$1]=$2}else{if($1==chr && $3==chr){d=int(sqrt(($2-$4)*($2-$4))*res); print $1,$2-o,$3,$4-o,$5/e[d],$6};}}' ${indecay} <( grep -v "#" ${inmatrix} ) >> ${outmatrix}
    done
done
