
scriptsdir=../../Analysis/
inputHiC=${scriptsdir}/cumulative_normalized_HiC_matrices_A_thaliana/

rep=$1
d=$2

wdir=contact_map    
outfilecorr=${PWD}/correlations_of_genomewide_contact_maps_in_genomewide_models_replica_${rep}_${d}.log


tag=Contact_map_genomewide
flag=intra
for res in 30kb ;
do
    resolution=$(echo $res | sed "s/kb/000/g")
    echo "Resolution ${resolution}"
    
    for dir in ${d} ;
    do
	echo $dir
	cd $dir
	
	version=$(echo $dir | sed "s/prodruns_//g")
	echo $version
	cd $wdir
	
	gzfile=$(ls -1 *gz | head -1)
	echo $gzfile
	if [[ -e ${gzfile} ]];
	then

	    outfile=${inputHiC}/cumulative_norm_HiC_matrices_full_at_${res}.tab
	    echo $outfile
	    awk '{print $2,$4,$5,$6}' <( awk '{if(NF!=0) print $0}' ${outfile} | grep -v "#") | awk '{if(NR==1){o1=$1;o2=$2};print $1-o1,$2-o2,$3}' > _exp_${rep}
	    wc -l _exp_${rep}
	    
	    for rc in 6.6667 ;
	    do
		if [[ ${rep} -eq 0 ]];
		then
		    #check=$(awk -v v=${version} -v t=${tag} -v rc=${rc} -v res=${resolution} 'BEGIN{c=0}{if($1==v && $2==t && $3==rc && $4==res && $6=="all" $NF==0){c++}}END{print c}' ${outfilecorr} )
		    check=0
		    if [[ ${check} -eq 0 ]]
		    then
			outfile2=contact_matrix_a_thaliana_at_${rc}sigma_at_${res}_${version}.tab
			if [[ -e ${outfile2} ]];
			then			    
			    echo $outfile2
			    grep -v "#" ${outfile2} | awk '{ print $2,$4,$5}' | awk '{if(NR==1){o1=$1;o2=$2};print $1-o1,$2-o2,$3}' > _model_${rep}
			    wc -l _model_${rep}

			    size=$(tail -1 _model_${rep} | awk '{if($1>$2){print $1+1}else{print $2+1}}')
			    echo $size
			    dcutoff=$(awk -v rc=$rc 'BEGIN{print int(rc*30)}')
			    corr=$(python ${scriptsdir}/compute_correlation_between_models_and_HiC_matrices.py $rep $size $dcutoff $flag) # 2> /dev/null)
			    
			    echo ${version} ${tag} ${rc} ${resolution} ${corr} "all" ${rep} >> ${outfilecorr}
			    rm _model_${rep}
			fi
		    fi
		fi # Condition on replica -eq 0

		if [[ ${rep} -gt 0 ]];
		then
		    for replica in ${rep} ;
		    do
			check=$(awk -v r=${rep} -v v=${version} -v t=${tag} -v rc=${rc} -v res=${resolution} 'BEGIN{c=0}{if($1==v && $2==t && $3==rc && $4==res && $6=="all" $NF==r){c++}}END{print c}' ${outfilecorr} )
			if [[ ${check} -eq 0 ]]
			then
			    mkdir  -p _tmp_replica_${replica}
			    cd _tmp_replica_${replica}
			    echo "Replica ${replica} $check"
			    outfile2=../contact_matrix_a_thaliana_at_${rc}sigma_at_${res}_replica_${replica}_${version}.tab
			    if [[ -e ${outfile2} ]];
			    then
				cp ../_exp_${replica} .
				grep -v "#" ${outfile2} | awk '{print $2,$4,$5}' > _model_${rep}
				size=$(tail -1 _model_${rep} | awk '{if($1>$2){print $1+1}else{print $2+1}}')
				dcutoff=$(awk -v rc=$rc 'BEGIN{print int(rc*30)}')
				corr=$(python ${scriptsdir}/compute_correlation_between_models_and_HiC_matrices.py $replica $size $dcutoff $flag 2> /dev/null)
				echo ${version} ${tag} ${rc} ${resolution} ${corr} "all" ${replica} >> ${outfilecorr}
			    fi
			    cd ..
			    rm -fr _tmp_replica_${replica}
			fi
		    done # Close cycle over $replica
		fi
	    done # Close cycle over $rc
	    rm -fr _exp_${rep}
	fi
	
	cd .. # Exit ${wdir}
	cd .. # Exit ${dir}
    done # Close cycle over $dir
done # Close cycle over $res
rm -fr _exp_${rep}

#for replica in $(seq 1 1 50);
#do
#done
