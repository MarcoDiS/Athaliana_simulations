scriptsdir=../../Analysis/
input=${scriptsdir}/contact_maps/

wdir=contact_map
indir=input_files

maxrc=6.6667
maxrct=6_6667
resolution=3000
format=tab
strcluster=IHIs

tmpIHIs=_tmp_IHIs
tmpfile=_tmp_add_${strcluster}
tmpstatefile=_tmp_states_${strcluster}
nIHIs=10
rm -fvr ${tmpstatefile} ${tmpfile}

#nreplicas=10

for dir in XXXwdirXXX ;
do
    echo $dir
    cd $dir #2> /dev/null
    pwd

    version=$(echo $dir | sed "s/prodruns_//g")
    echo $version

    mkdir -p ${wdir}
    cd ${wdir}

    echo "Calculating the contacts within ${maxrc}sigma distance"
    if [[ -e distances_model_a_thaliana_replica_1_below_${maxrct}sigma_at_3kb_txt.tar.gz ]];
    then
	for rc in 6.6667 ;
	do

	    echo "dcutoff = ${rc}"
	    if [[ ! -e ${tmpIHIs} ]];
	    then
		cat <( zcat distances_model_a_thaliana_replica_*_below_${maxrct}sigma_at_3kb_txt.tar.gz ) | awk -v rc=${rc} -v s=${size} '{if($4<=(rc*rc) && $5>=1000000 && $5<=20000000){print $0}}' > ${tmpIHIs}
	    fi

	    for res in 3kb ;
	    do		
		for nl1 in XXXnIHI1XXX; #$(seq 1 1 ${nIHIs});
		do
		    IHI1=$(awk -v nl1=$nl1 '{if(NR==nl1) print $1}' <(cat ${scriptsdir}/Utilities/IHIs_positions.txt | grep A | grep -v "#")) ;
		    cbead1=$(grep -w ${IHI1} ${scriptsdir}/Utilities/IHIs_positions.txt | awk '{print int(($2+$3)/2)}')
		    centralbead1=$(awk -v cb=${cbead1} '{if($1==cb)print $2}' ${scriptsdir}/Utilities/bin_mapping_at_3kb.txt)
		    awk -v cb=${cbead1} 'BEGIN{print cb}{if($1==cb)print $2}' ${scriptsdir}/Utilities/bin_mapping_at_3kb.txt
		    
		    # 166beads = 498000bp
		    startbead1=$(awk -v cb=${centralbead1} 'BEGIN{print cb-166}')
		    stopbead1=$(awk -v cb=${centralbead1} 'BEGIN{print cb+166}')
		    echo "ID: $IHI1 ; cbead: $cbead1 ; startbead: $startbead1 ; centralbead: $centralbead1 ; stopbead: ${stopbead1}"

		    for nl2 in XXXnIHI2XXX; #$(seq $((${nl1}+1)) 1 ${nIHIs});
		    do
			IHI2=$(awk -v nl2=$nl2 '{if(NR==nl2) print $1}' <(cat ${scriptsdir}/Utilities/IHIs_positions.txt | grep A | grep -v "#")) ;
			cbead2=$(grep -w ${IHI2} ${scriptsdir}/Utilities/IHIs_positions.txt | awk '{print int(($2+$3)/2)}')
			centralbead2=$(awk -v cb=${cbead2} '{if($1==cb)print $2}' ${scriptsdir}/Utilities/bin_mapping_at_3kb.txt)

    		        # 166beads = 498000bp
			startbead2=$(awk -v cb=${centralbead2} 'BEGIN{print cb-166}')
			stopbead2=$(awk -v cb=${centralbead2} 'BEGIN{print cb+166}')
			echo "ID: $IHI2 ; cbead: $cbead2 ; startbead: $startbead2 ; centralbead: $centralbead2 ; stopbead: ${stopbead2}"

			outfilepair=contact_matrix_a_thaliana_at_${rc}sigma_at_${res}_${version}_around_${IHI1%A}_and_${IHI2%A}.${format}
			if [[ ! -s ${outfilepair} ]];
			then
			    echo $IHI1 $IHI2
			    
			    size=$((${stopbead1}-${startbead1}+1))
			    echo $size
			    
			    awk -v s=$size -v startbead1=$startbead1 -v stopbead1=$stopbead1 -v startbead2=$startbead2 -v stopbead2=$stopbead2 'BEGIN{for(i=0;i<s;i++)for(j=0;j<s;j++)m[i,j]=0}{if((startbead1<=$2 && $2<=stopbead1) && (startbead2<=$3 && $3<=stopbead2)){ind1=$2-startbead1; ind2=$3-startbead2; m[ind1,ind2]+=1}; if((startbead1<=$3 && $3<=stopbead1) && (startbead2<=$2 && $2<=stopbead2)){ind1=$3-startbead1; ind2=$2-startbead2; m[ind1,ind2]+=1}}END{for(i=0;i<s;i++)for(j=0;j<s;j++){print i,j,m[int(i),int(j)]}}' ${tmpIHIs} > ${outfilepair}
			fi
		    done # Close cycle over $IHI2
		done # Close cycle over $IHI1
		

		outfile=contact_matrix_a_thaliana_at_${rc}sigma_at_${res}_${version}_around_IHIs.${format}
		if [[ ! -s ${outfile} ]];
		then
		    echo "Size: $size ; Outfile: ${outfile}"
		    
		    echo "# CRM 1Mb around IHIs regions centred at the bin 166" > ${outfile}
		    echo "#parti partj Obs " >> ${outfile}
			
		    awk -v s=${size} '{m[int($1),int($2)]+=$3}END{for(i=0;i<s;i++)for(j=0;j<s;j++)print i,j,m[i,j]/45}' contact_matrix_a_thaliana_at_${rc}sigma_at_${res}_${version}_around_IHI*_and_IHI*.${format} >> ${outfile}

		fi
		    
		tail -3 ${outfile} 
		
	    done
	    #rm -fvr _tmp
	done
    fi
    
    cd .. # Exit ${wdir}
    cd .. # Exit ${dir}

done
echo "Analysis DONE!"
