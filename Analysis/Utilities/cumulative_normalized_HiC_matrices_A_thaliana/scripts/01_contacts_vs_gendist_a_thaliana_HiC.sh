# Compute contacts vs L
for res in 30kb ;
do
    echo $res
    resolution=$(echo ${res} | sed 's/kb/000/g')
    infile=cumulative_norm_HiC_matrices_full_at_${res}.tab

    # Genome-wide
    echo "#Gendist avg_contacts for genome-wide matrix"
    outfile=contacts_vs_gendist_a_thaliana_norm_full_${res}.txt

    awk -v res=${resolution} '{if($1==$3){d=sqrt(($2-$4)*($2-$4));if(d!=0) h[d]+=$5; h2[d]+=$5*$5; cnt[d]++};}END{factor=h[10]/cnt[10]; for(i in h){avg=h[i]/(cnt[i]*factor); avg2=h2[i]/(cnt[i]*factor*factor); stddev=sqrt(avg2-avg*avg); print (i)*res,avg,stddev}}' <( grep -v MASKED ${infile} ) | sort -k 1n > ${outfile}

    # Per chromosome    
    for chr in chr1 chr2 chr3 chr4 chr5 ;
    do
	c=$(echo $chr | sed "s/chr//g" | awk '{print $1"_dna"}')	
	outfile=contacts_vs_gendist_a_thaliana_norm_${c}_at_${res}.txt

	echo "#Gendist avg_contacts for ${chr}"
	    awk -v res=${resolution} -v chr=${chr} '{if($1==chr && $3==chr){d=sqrt(($2-$4)*($2-$4));if(d!=0) h[d]+=$5; h2[d]+=$5*$5; cnt[d]++};}END{factor=h[10]/cnt[10]; for(i in h){avg=h[i]/(cnt[i]*factor); avg2=h2[i]/(cnt[i]*factor*factor); stddev=sqrt(avg2-avg*avg); print (i)*res,avg,stddev}}' <( grep -v MASKED ${infile} ) | sort -k 1n > ${outfile}

    done
done