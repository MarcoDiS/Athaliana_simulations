
scriptsdir=../../Analysis/
inputHiC=${scriptsdir}/Utilities/cumulative_normalized_HiC_matrices_A_thaliana/

rep=$1
d=$2

wdir=contact_map    
outfilecorr=${PWD}/correlations_of_intra_chromosome_contact_maps_in_genomewide_models_replica_${rep}_${d}.log


for chr1 in 1 2 3 4 5 ;
do
    echo "Chr ${chr1}"

    for chr2 in 1 2 3 4 5 ;
    do
	echo "Chr ${chr2}"
	c1=$(echo ${chr1} | sed "s/chr//g")
	c2=$(echo ${chr2} | sed "s/chr//g")
	if [[ ${c2} -gt ${c1} ]];
	then
	    continue
	fi 
	tag="Contact_map_inter"
	flag="inter"
	if [[ ${c1} -eq ${c2} ]];
	then
	    tag="Contact_map_intra"
	    flag="intra"
	#else
	    #continue
	fi
	

	for res in 30kb ;
	do
	    resolution=$(echo $res | sed "s/kb/000/g")
	    echo "Resolution ${resolution}"
	    
	    for dir in ${d} ;
	    do

		outfile=${inputHiC}/cumulative_norm_HiC_matrices_full_at_${res}.tab
		echo $outfile
		echo "Exp matrix"
		awk -v c1="chr"${chr1} -v c2="chr"${chr2} '{if($1==c1 && $3==c2) print $2,$4,$5,$6}' <( awk '{if(NF!=0) print $0}' ${outfile} | grep -v "#") | awk '{if(NR==1){o1=$1;o2=$2};print $1-o1,$2-o2,$3}' > ${dir}/${wdir}/_exp_${rep}
		awk -v c1="chr"${chr1} -v c2="chr"${chr2} '{if($1==c1 && $3==c2) print $2,$4,$5,$6}' <( awk '{if(NF!=0) print $0}' ${outfile} | grep -v "#") | awk '{if(NR==1){o1=$1;o2=$2; print o1,o2}}'
		sort -k 1n _exp_${rep} | uniq | head
		sort -k 1n _exp_${rep} | uniq | tail
		sort -k 2n _exp_${rep} | uniq | head
		sort -k 2n _exp_${rep} | uniq | tail
		
		echo $dir
		cd $dir
		
		version=$(echo $dir | sed "s/prodruns_//g")
		echo $version
		cd $wdir
		
		#cp ../../_exp_${rep} .
		
		for rc in 6.6667 ; #1.5 3.3334 6.6667 ;
		do
		    if [[ ${rep} -eq 0 ]];
		    then
			#cp _exp_${rep} _exp_${rep}
			check=$(awk -v v=${version} -v t=${tag} -v rc=${rc} -v res=${resolution} -v chr1=${chr1} -v chr2=${chr2} 'BEGIN{c=0}{if($1==v && $2==t && $3==rc && $4==res && $6==chr1"-"chr2 && $NF==0){c++}}END{print c}' ${outfilecorr} )
			if [[ ${check} -eq 0 ]]
			then
			    outfile2=contact_matrix_a_thaliana_at_${rc}sigma_at_${res}_${version}.tab
			    if [[ -e ${outfile2} ]];
			    then			    
				echo $outfile2
				echo "Model matrix"
				grep -v "#" ${outfile2} | awk -v chr1=${chr1} -v chr2=${chr2} '{if($1=="chr"chr1 && $3=="chr"chr2) print $2,$4,$5}' | awk '{if(NR==1){o1=$1;o2=$2};print $1-o1,$2-o2,$3}' > _model_${rep}
				#grep -v "#" ${outfile2} | awk -v chr1=${chr1} -v chr2=${chr2} '{if($1=="chr"chr1 && $3=="chr"chr2) print $2,$4,$5}' | awk '{if(NR==1){o1=$1;o2=$2; print o1,o2}; if($1-o1<0 || $2-o2<0){print $0,NR; exit}}'
				sort -k 1n _model_${rep} | uniq | head 
				sort -k 1n _model_${rep} | uniq | tail 
				sort -k 2n _model_${rep} | uniq | head 
				sort -k 2n _model_${rep} | uniq | tail

				wc -l _model_${rep}				
				size=$(tail -1 _model_${rep} | awk '{if($1>$2){print $1+1}else{print $2+1}}')
				echo "Size: $size"
				dcutoff=$(awk -v rc=$rc 'BEGIN{print int(rc*30)}')
				#replica  = int(sys.argv[1])
				#size    = int(sys.argv[2])
				#dcutoff = float(sys.argv[3])
				#flag    = sys.argv[4]
				head _model_${rep} _exp_${rep}
				corr=$(python ${scriptsdir}/compute_correlation_between_models_and_HiC_matrices.py $rep $size $dcutoff $flag)
				
				echo ${version} ${tag} ${rc} ${resolution} ${corr} "${chr1}-${chr2}" ${rep} >> ${outfilecorr}
				    #mv plot.png correlation_plots_chr${chr1}_chr${chr2}_at_${rc}sigma_at_${res}_${version}.png
				rm _model_${rep} plot.pdf
			    fi
			fi
		    fi # Condition on replica -eq 0

		    if [[ ${rep} -gt 0 ]];
		    then
			for replica in ${rep} ;
			do
			    check=$(awk -v v=${version} -v t=${tag} -v rc=${rc} -v rep=${replica} -v res=${resolution} -v chr1=${chr1} -v chr2=${chr2} 'BEGIN{c=0}{if($1==v && $2==t && $3==rc && $4==res && $6==chr1"-"chr2 && $NF==rep){c++}}END{print c}' ${outfilecorr} )	
			    if [[ ${check} -eq 0 ]]
			    then
				mkdir  -p _tmp_replica_${replica}
				cd _tmp_replica_${replica}
				echo "Replica ${replica} $check"
				outfile2=../contact_matrix_a_thaliana_at_${rc}sigma_at_${res}_replica_${replica}_${version}.tab
				if [[ -e ${outfile2} ]];
				then
				    cp ../_exp_${rep} .
				    echo "Model matrix"
				    grep -v "#" ${outfile2} | awk -v chr1=${chr1} -v chr2=${chr2} '{if($1=="chr"chr1 && $3=="chr"chr2) print $2,$4,$5}' | awk '{if(NR==1){o1=$1;o2=$2};print $1-o1,$2-o2,$3}' > _model_${rep}
				    sort -k 1n _model_${rep} | uniq
				    sort -k 2n _model_${rep} | uniq				    

				    size=$(tail -1 _model_${rep} | awk '{if($1>$2){print $1+1}else{print $2+1}}')
				    echo "Size: $size"
				    dcutoff=$(awk -v rc=$rc 'BEGIN{print int(rc*30)}')
				    corr=$(python ${scriptsdir}/compute_correlation_between_models_and_HiC_matrices.py $rep $size $dcutoff $flag)
				    echo ${version} ${tag} ${rc} ${resolution} ${corr} "${chr1}-${chr2}" ${rep} >> ${outfilecorr}
				fi
				cd ..
				rm -fr _tmp_replica_${replica}
			    fi
			done # Close cycle over $replica
		    fi
		done # Close cycle over $rc
		rm -fr _exp_${rep}
		
		cd .. # Exit ${wdir}
		cd .. # Exit ${dir}
	    done # Close cycle over $dir
	done # Close cycle over $res
    done # Close cycle over $chr2
done # Close cycle over $chr1
rm -fr _exp_${rep}

