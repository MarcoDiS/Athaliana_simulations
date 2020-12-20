scriptsdir=../../Analysis/

inputHiC=${scriptsdir}/Utilities/cumulative_normalized_HiC_matrices_A_thaliana/

rep=$1
d=$2

wdir=contact_map    
outfilecorr=${PWD}/correlations_of_genomewide_contacts_vs_L_in_genome_wide_models_replica_${rep}_${d}.log

tag=Contacts_vs_L_genomewide

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
	
	for rc in 6.6667 ;
	do

	    if [[ -d  _tmp_replica_${rep} ]];
	    then
		exit
	    fi

	    mkdir  -p _tmp_replica_${rep}
	    cd _tmp_replica_${rep}

            touch ${outfilecorr}

	    check=$(awk -v r=${rep} -v v=${version} -v t=${tag} -v rc=${rc} -v res=${resolution} 'BEGIN{c=0}{if($1==v && $2==t && $3==rc && $4==res && $7==r){c++}}END{print c}' ${outfilecorr} 2> /dev/null)
	    if [[ ${check} -ne 0 ]]
	    then
		echo "Calculation for replica ${rep} for ${version} of ${tag} at ${rc}sigma ${resolution}bp for genome-wide DONE!"
		cd .. # Exit from _tmp_replica_${rep} 
		continue
	    fi

	    outfile=${inputHiC}/contacts_vs_gendist_a_thaliana_norm_full_${res}.txt
	    echo $outfile
	    awk '{print $2}' ${outfile} > _exp_${rep}
	    wc -l _exp_${rep}

	    if [[ ${rep} -eq 0 ]];
	    then
		outfile2=../contacts_vs_gendist_a_thaliana_at_${rc}sigma_at_${res}_${version}.txt
		if [[ -e ${outfile2} ]];
		then			    
		    echo $outfile2
		    awk '{if($1>0) print $2}' ${outfile2} > _model_${rep}
		    wc -l _model_${rep}
		    
		    size=$(tail -1 _model_${rep} | awk '{if($1>$2){print $1+1}else{print $2+1}}')
		    echo $size
		    echo "Computing the correlation"
		    dcutoff=$(awk -v rc=$rc 'BEGIN{print int(rc*30)}')
		    correl=$(python ${scriptsdir}/Utilities/compute_correlation_between_HiC_and_models_contact_maps.py ${rep})

		    echo $correl

		    corr=$(echo ${correl} | awk '{print $1}')
                    npoints=$(echo ${correl} | awk '{print $2}')
		    
		    echo ${version} ${tag} ${rc} ${resolution} ${corr} "all" ${rep} ${npoints} >> ${outfilecorr}
		fi

	    fi # Condition on replica -eq 0
	    
	    if [[ ${rep} -gt 0 ]];
	    then
		for replica in ${rep} ;
		do
		    echo "Replica ${replica} $check"
		    outfile2=../contacts_vs_gendist_a_thaliana_at_${rc}sigma_at_${res}_replica_${replica}_${version}.txt
		    if [[ -e ${outfile2} ]];
		    then
			echo $outfile2
			awk '{if($1>0) print $2}' ${outfile2} > _model_${replica}
			wc -l _model_${rep}
			
			size=$(tail -1 _model_${replica} | awk '{if($1>$2){print $1+1}else{print $2+1}}')
			echo $size
			
			echo "Computing the correlation"
			dcutoff=$(awk -v rc=$rc 'BEGIN{print int(rc*30)}')
			correl=$(python ${scriptsdir}/Utilities/compute_correlation_between_HiC_and_models_contact_maps.py ${rep})
		
			corr=$(echo ${correl} | awk '{print $1}')
                        npoints=$(echo ${correl} | awk '{print $2}')
	
			echo $correl

			echo ${version} ${tag} ${rc} ${resolution} ${corr} "all" ${replica} ${npoints} >> ${outfilecorr}
		    fi
		done # Close cycle over $replica
	    fi
	    cd .. # Exit _tmp_replica_${replica}
	done # Close cycle over $rc
	rm -fvr _tmp_replica_${rep}
	
	cd .. # Exit ${wdir}
	cd .. # Exit ${dir}
    done # Close cycle over $dir
done # Close cycle over $res
