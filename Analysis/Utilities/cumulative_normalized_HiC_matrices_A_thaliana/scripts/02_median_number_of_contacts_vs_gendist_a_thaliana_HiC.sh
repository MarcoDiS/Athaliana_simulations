# Compute contacts vs L
for res in 30kb ;
do
    echo $res
    resolution=$(echo ${res} | sed 's/kb/000/g')
    infile=cumulative_norm_HiC_matrices_full_at_${res}.tab

    # Genome-wide
    (
	echo "#Gendist avg_contacts for genome-wide matrix"
	outfile=median_number_of_contacts_vs_gendist_a_thaliana_norm_full_${res}.txt
	
	rm -fvr ${outfile}
	tmpfile=_tmp_sorted_full
	awk -v res=${resolution} '{if($1==$3){d=sqrt(($2-$4)*($2-$4)); print d*res,$5}else{print "Inter",$5}}' <( grep -v MASKED ${infile} ) | sort -k 2n,2n > ${tmpfile}
	for d in $(awk '{print $1}' ${tmpfile} | uniq | sort -k 1n | uniq);
	do
	    awk -v d=$d 'BEGIN{cnt=0}{if($1==d){cnt++;v[cnt]=$2}}END{if(cnt%2==1){print d,v[int(cnt*0.5)+1]}else{print d,(v[int(cnt*0.5)+1]+v[int(cnt*0.5)])*0.5}}' ${tmpfile} >> ${outfile}
	    tail -1 ${outfile}
	done
    ) &

    # Per chromosome
    for chr in chr1 chr2 chr3 chr4 chr5 ;
    do
	(
	    c=$(echo $chr | sed "s/chr//g" | awk '{print $1"_dna"}')

	    echo "#Gendist avg_contacts for ${chr}"	    
	    outfile=median_number_of_contacts_vs_gendist_a_thaliana_norm_${c}_at_${res}.txt

	    rm -fvr ${outfile}
	    tmpfile=_tmp_sorted_${chr}
	    awk -v res=${resolution} -v chr=${chr} '{if($1==chr && $1==$3){d=sqrt(($2-$4)*($2-$4)); print d*res,$5}else{print "Inter",$5}}' <( grep -v MASKED ${infile} ) | sort -k 2n,2n > ${tmpfile}
	    for d in $(awk '{print $1}' ${tmpfile} | uniq | sort -k 1n | uniq);
	    do
		awk -v d=$d 'BEGIN{cnt=0}{if($1==d){cnt++;v[cnt]=$2}}END{if(cnt%2==1){print d,v[int(cnt*0.5)+1]}else{print d,(v[int(cnt*0.5)+1]+v[int(cnt*0.5)])*0.5}}' ${tmpfile} >> ${outfile}
		tail -1 ${outfile}
	    done
	) &
    done
done
