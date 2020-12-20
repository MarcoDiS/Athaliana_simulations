resolution=3000 # bp per bead
BinWidth=10 # Number of beads per bin
reskb=$(awk -v res=${resolution} -v bw=${BinWidth} 'BEGIN{print int(res*bw/1000)"kb"}')

infile=chromatin_states_haploid_genome_with_locations.txt
outfile=chromatin_states_haploid_genome_with_locations_at_${reskb}_new.txt
rm -fr ${outfile}

#head -3 chromatin_states_haploid_genome_with_locations.txt
#tail -3 chromatin_states_haploid_genome_with_locations.txt
stopbin=-1
for chr in chr1 chr2 chr3 chr4 chr5 ;
do
    cat ${infile} | grep ${chr} | grep -v NOR | head -3
    cat ${infile} | grep ${chr} | grep -v NOR | tail -3
    echo

    totbins=$(awk -v chr=${chr} 'BEGIN{min=100000}{if($1==chr && $3==chr){if($2>max) max=$2; if($2<min){min=$2}; if($4>max) max=$4; if($4<min){min=$4}}}END{print min" - "max}' /home/devel/mstefano/cumulative_normalized_HiC_matrices_A_thaliana/cumulative_norm_HiC_matrices_full_at_30kb.tab)
    echo "Bins in Hi-C maps for ${chr} = ${totbins}"

    nbins=$(cat ${infile} | grep -v NOR | grep ${chr} | wc -l)
    startbin=$((${stopbin}+1))
    stopbin=$(awk -v s1=${stopbin} -v nb=${nbins} -v bw=${BinWidth} 'BEGIN{print s1+int(nb/bw)+1}')
    #echo "Nbins = ${nbins}"
    if [[ ${chr} == "chr3" ]];
    then
	stopbin=$((${stopbin}-1))
    fi
    echo "Bins in state assignment for ${chr} = ${startbin} - ${stopbin}"
    echo

    for bin in $(seq ${startbin} 1 ${stopbin});
    do
	start=$(awk -v sb=${startbin} -v bw=${BinWidth} -v b=${bin} -v res=${resolution} 'BEGIN{print (b-sb)*res*bw}')
	stop=$(awk -v sb=${startbin} -v bw=${BinWidth} -v b=${bin} -v res=${resolution} 'BEGIN{print (b-sb+1)*bw*res}')
	
	state=$(sed -e "s/:/ /g" -e "s/-/ /g" ${infile} | grep -v NOR | grep ${chr} | awk -v s1=${start} -v s2=${stop} 'BEGIN{flag=0}{if($2==s1){flag=1}; if(flag==1){print $0}; if($3==s2){flag=0}}' | awk '{h[$NF]++}END{for(i in h) print i,h[i]}' | sort -g -k2n | tail -1 | awk '{print $1}')

	mask=$(sed -e "s/:/ /g" -e "s/-/ /g" ${infile} | grep -v NOR | grep ${chr} | awk -v s1=${start} -v s2=${stop} 'BEGIN{flag=0}{if($2==s1){flag=1}; if(flag==1){print $0}; if($3==s2){flag=0}}' | awk '{h[$NF]++}END{for(i in h) print i,h[i]}' | sort -g -k2n | tac | awk '{h[$2]++; if($2>=max){max=$2;if(states==""){states=$1}else{states=states"_"$1}}}END{if(h[max]>1){print "NOT-ASSIGNED"}else{print "OK"}}') 

	if [[ $mask == "NOT-ASSIGNED" ]];
	then
	    state="not-assigned"
	fi

	flag=$(sed -e "s/:/ /g" -e "s/-/ /g" ${infile} | grep -v NOR | grep ${chr} | awk -v s1=${start} -v s2=${stop} 'BEGIN{flag=0}{if($2==s1){flag=1}; if(flag==1){print $0}; if($3==s2){flag=0}}' | awk '{h[$NF]++}END{for(i in h) print i,h[i]}' | sort -g -k2n | tac | awk '{h[$2]++; if($2>=max){max=$2;if(states==""){states=$1}else{states=states"_"$1}}}END{if(h[max]>1){print "CHECK_"states}else{print "OK"}}') 

	sed -e "s/:/ /g" -e "s/-/ /g" ${infile} | grep -v NOR | grep ${chr} | awk -v s1=${start} -v s2=${stop} 'BEGIN{flag=0}{if($2==s1){flag=1}; if(flag==1){print $0}; if($3==s2){flag=0}}' | awk -v b=${bin} -v f=${flag} '{h[$NF]++}END{for(i in h) print "#"f,b,i,h[i]}' | sort -g -k4n	>> ${outfile}
	
	echo "${chr}:${start}-${stop} ${bin} ${state} ${mask}"
	echo "${chr}:${start}-${stop} ${bin} ${state}" >> ${outfile}

    done
    cat ${outfile} | grep ${chr} | head -3
    cat ${outfile} | grep ${chr} | tail -3
    echo
    
done