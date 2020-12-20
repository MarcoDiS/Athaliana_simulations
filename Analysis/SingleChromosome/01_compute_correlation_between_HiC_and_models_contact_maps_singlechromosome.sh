scriptsdir=../../Analysis/
input=${scriptsdir}/contact_maps/

wdir=contact_map
indir=input_files
format=tab

nreplicas=10

for dir in XXXsystemXXX ;
do
    echo $dir
    cd $dir 2> /dev/null

    version=$(echo $dir | sed "s/prodruns_//g")

    if [[ -d ${wdir} ]];
    then
	cd $wdir
	ls -lrth

	for replica in $(seq 1 1 ${nreplicas})
	do
	    
	    gzfile=distances_model_a_thaliana_replica_${replica}_below_6_6667sigma_at_3kb_txt.tar.gz
	    ls -lrth $gzfile

	    if [[ -e ${gzfile} ]];
	    then
		if [[ $( ls -1 _tmp_replica_${replica} 2> /dev/null | wc -l | awk '{print $1}') -gt 0 ]]; 
		then 
		    exit
		fi
		
		cat <( zcat ${gzfile} ) | awk '{if(NF==5){f="_tmp_replica_"$1; print $0 > f}}'		
		
		for rc in 6.6667 ;
		do
		    
		    ntab=$(ls -1 *${rc}*replica_${replica}_*tab 2> /dev/null | wc -l | awk '{print $1}')
		    if [[ ${ntab} -ge 3 ]];
		    then
			continue
		    fi

		    echo "dcutoff = ${rc}"
	    	    awk -v rc=${rc} '{if($4<=(rc*rc) && $NF>=1000000) print $0}' _tmp_replica_${replica} > _tmp_replica_${replica}_${rc}
		
		    for res in 30kb ; #150kb 30kb 15kb ;
		    do
                        # Produce genome-wide model contact map
			resolution=$(echo ${res} | sed "s/kb/000/g")
		        # Bins of the NOR to remove
			MappingFile=bin_mapping_at_${res}_chr4.txt		
			echo ${resolution}
			if [[ ! -e ${scriptsdir}/${MappingFile} ]];
			then
			    echo "ERROR! File ${scriptsdir}/${MappingFile} missing"
			    echo "Create it first using ${scriptsdir}/bin_mapping.sh"
			    continue
			fi
			size=$(tail -1 ${scriptsdir}/${MappingFile} | awk '{print $4+1}')		    

			outfile=contact_matrix_a_thaliana_at_${rc}sigma_at_${res}_replica_${replica}_${version}.${format}
			outfiledecay=contacts_vs_gendist_a_thaliana_at_${rc}sigma_at_${res}_replica_${replica}_${version}.txt
			outfilenormdecay=norm_contacts_vs_gendist_a_thaliana_at_${rc}sigma_at_${res}_replica_${replica}_${version}.txt
			echo $size
			
			if [[ ! -e ${outfile} ]];
			then
			    echo ${outfile}
			    
			    if [[ ${format} == "tab" ]];
			    then
				
				awk -v s=${size} 'BEGIN{for(i=0;i<s;i++)for(j=0;j<=i;j++)print 11,i,j,0,0}' > _tmp_replica_add
				
				echo "# CRM 4_dna 18585056" >> ${outfile}
				echo "#chri parti chrj partj Obs Obs/Exp" >> ${outfile}
				
				cat <(cat ${scriptsdir}/${MappingFile} _tmp_replica_${replica}_${rc} | awk '{if(NF==4){b[$1-1]=$4; c[$1-1]=$2}; if(NF==5){if($2 in b) if($3 in b){ind1=b[$2] ; chr1=c[$2] ; ind2=b[$3] ; chr2=c[$3] ; if(ind2<=ind1){print chr1"_"ind1"_"chr2"_"ind2}else{print chr2"_"ind2"_"chr1"_"ind1}}}}') <(cat ${scriptsdir}/${MappingFile} _tmp_replica_add | awk '{if(NF==4){b[$4]=$4; c[$4]=$2}; if(NF==5){if($2 in b) if($3 in b){ind1=b[$2] ; chr1=c[$2] ; ind2=b[$3] ; chr2=c[$3] ; if(ind2<=ind1){print chr1"_"ind1"_"chr2"_"ind2}else{print chr2"_"ind2"_"chr1"_"ind1}}}}') | awk '{h[$1]++}END{for(i in h) print i,h[i]-1}' | sed "s/_/ /g" | sort -k 2n,2n -k 4n,4n > _tmp_replica_matrix
				
				awk '{if($1==$3){d=$2-$4;print d,$5}}' _tmp_replica_matrix | awk '{e[$1]+=$2;e2[$1]+=($2*$2);c[$1]++}END{for(i in e){avg=e[i]/c[i];avg2=e2[i]/c[i];stddev=sqrt(avg2-avg*avg);print i,avg,stddev}}' > _tmp_replica_decay
				cat _tmp_replica_decay _tmp_replica_matrix | awk '{if(NF==3){w[$1]=$2};if(NF==5){corrected="NA";if($1==$3){d=$2-$4; if(w[d]!=0){corrected=$5/w[d]}else{corrected=$5}}; print $0,corrected}}' >> ${outfile}
				
				awk -v r=${resolution} '{c[$1*r]=$2;e[$1*r]=$3;f=1}END{for(i in e) print i,c[i]/f,e[i]/f}' _tmp_replica_decay | sort -k 1n > ${outfiledecay}
				awk -v r=${resolution} '{c[$1*r]=$2;e[$1*r]=$3;if($1*r==300000)f=$2}END{for(i in e) print i,c[i]/f,e[i]/f}' _tmp_replica_decay | sort -k 1n > ${outfilenormdecay}
				
				rm -fvr _tmp_replica_matrix _tmp_replica_decay _tmp_replica_add
			    fi
			    
			    tail -3 ${outfile} #| awk '{print $(NF-2),$(NF-1),$NF}'
			    
			fi

		    done # Close cycle over $res
		    rm _tmp_replica_${replica}_${rc}
		done  # Close cycle over $rc
		rm _tmp_replica_${replica}
	    fi # Close conditional over ${gzfile}
	    
	done  # Close cycle over $replica
	rm -fvr _tmp_replica_${replica}
	
	cd .. # Exit ${wdir}
    fi # Close condistional over ${wdir}
    cd .. # Exit ${dir}
    
done # Close cycle over $dir
echo "Analysis DONE!"
