scriptsdir=../../Analysis/
input=${scriptsdir}/contact_maps/

wdir=contact_map
indir=input_files

maxrc=6.6667
maxrct=6_6667
resolution=3000

for dir in XXXwdirXXX ;
do
    echo $dir
    cd $dir #2> /dev/null
    pwd

    version=$(echo $dir | sed "s/prodruns_//g")
    echo $version

    mkdir -p ${wdir}
    cd ${wdir}
    tmpH3K27me3=_tmp_H3K27me3

    echo "Calculating the contacts within ${maxrc}sigma distance"
    if [[ -e distances_model_a_thaliana_replica_1_below_${maxrct}sigma_at_3kb_txt.tar.gz ]];
    then
	for rc in 6.6667 ;
	do

	    echo "dcutoff = ${rc}"
	    if [[ ! -e ${tmpH3K27me3} ]];
	    then
		cat <( zcat distances_model_a_thaliana_replica_*_below_${maxrct}sigma_at_3kb_txt.tar.gz ) | awk -v rc=${rc} -v s=${size} '{if($4<=(rc*rc) && $5>=1000000 && $5<=20000000){print $0}}' > ${tmpH3K27me3}
	    fi

	    for res in 3kb ;
	    do
		for cluster in XXXclusterXXX ;
		do		    
		    chr=$(grep Plot ${scriptsdir}/Utilities/H3K27me3_clusters.txt | grep -w ${cluster} | awk '{print $NF}' | sed "s/:/ /g" | awk '{print $1}')
		    if [[ ! -e _${chr} ]];
		    then
			grep ${chr} ${scriptsdir}/Utilities/chromatin_states_haploid_genome.txt | awk '{print $1,$2,$3}' > _${chr}
		    fi
		    #cat _$chr
		    start=$(grep Plot ${scriptsdir}/Utilities/H3K27me3_clusters.txt | grep -w ${cluster} | awk '{print $NF}' | sed -e "s/:/ /g" -e "s/\-/ /g" | awk '{print (int($2/3000))*3000}')
		    stop=$(grep Plot ${scriptsdir}/Utilities/H3K27me3_clusters.txt | grep -w ${cluster} | awk '{print $NF}' | sed -e "s/:/ /g" -e "s/\-/ /g" | awk '{print (int($3/3000))*3000}')
		    
		    startbead=$(awk -v res=$resolution -v start=${start} -v c=${chr} 'BEGIN{s=int(start/res)+1}{if($3!="NORs"){c++;if(c==s){print $2}}}' _${chr})
		    stopbead=$(awk -v res=$resolution -v stop=${stop} -v c=${chr} 'BEGIN{s=int(stop/res)}{if($3!="NORs"){c++;if(c==s){print $2}}}' _${chr})
		    echo $cluster $chr $start $startbead $stop ${stopbead}

		    size=$((${stopbead}-${startbead}+1))

		    format=tab
		    strcluster=$(echo $cluster | sed "s/:/_/g")
		    outfile=contact_matrix_a_thaliana_at_${rc}sigma_at_${res}_${version}_${strcluster}.${format}
		    tmpfile=_tmp_add_${strcluster}
		    tmpstatefile=_tmp_states_${strcluster}
		    tmpmatrixfile=_tmp_matrix_${strcluster}
		    #rm -fr ${outfile}
		    echo $size

		    #Check epigenetic state of the bead in the cluster
		    #chr5 42251 Telomere
		    grep ${chr} ${scriptsdir}/Utilities/chromatin_states_haploid_genome.txt | grep -v NORs | awk -v s1=$start -v s2=$stop -v res=$resolution '{if(int(s1/res)<NR && NR<=int(s2/res)) print $3}' > ${tmpstatefile}

		    cat ${tmpstatefile}
		    if [[ ! -e ${outfile} ]];
		    then
			echo ${outfile}
			
			if [[ ${format} == "tab" ]];
			then
			    
			    awk -v s=${size} '{if(NF==1){n++;state[n-1]=$1}}END{for(i=0;i<s;i++)for(j=0;j<=i;j++)print 51,i,state[i],j,state[j],0,0}' ${tmpstatefile} > ${tmpfile}
			    cat ${tmpfile}

			    echo "# CRM ${cluster} ${chr}:${start}-${stop}" >> ${outfile}
			    echo "#chri statei parti chrj statej partj Obs Obs/Exp" >> ${outfile}

			    awk -v s=${size} -v c=${chr} -v startbead=$startbead -v stopbead=$stopbead '{if($1==51){m[int($2),int($4)]=$6;state[int($2)]=$3;state[int($4)]=$5}else{if(((startbead-1)<=$2 && $2<=(stopbead-1)) && ((startbead-1)<=$3 && $3<=(stopbead-1))){ind1=$2; ind2=$3; if($3>$2){ind1=$3; ind2=$2}; m[int(ind1-(startbead-1)),int(ind2-(startbead-1))]+=1}}}END{for(i=0;i<s;i++){for(j=0;j<=i;j++){d=sqrt((i-j)*(i-j));e[d]+=m[i,j];cnt[d]++}}; for(i=0;i<s;i++)for(j=0;j<=i;j++){d=sqrt((i-j)*(i-j));print c,state[i],i,c,state[j],j,m[i,j],m[i,j]/((e[d]+1)/(cnt[d]+1))}}' ${tmpfile} ${tmpH3K27me3} >> ${outfile} 
			fi
			
			tail -3 ${outfile} #| awk '{print $(NF-2),$(NF-1),$NF}'
			rm ${tmpfile} ${tmpstatefile} ${tmpmatrixfile}
		    fi
		done # Close cycle over $cluster
	    done
	    #rm -fvr _tmp
	done
    fi
    
    cd .. # Exit ${wdir}
    cd .. # Exit ${dir}

done
echo "Analysis DONE!"
