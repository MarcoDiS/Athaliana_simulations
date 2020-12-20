scriptsdir=../../Analysis/
input=${scriptsdir}/Utilities/contact_maps/
inputHiC=${scriptsdir}/Utilities/cumulative_normalized_HiC_matrices_A_thaliana/

wdir=contact_map
indir=input_files

nreplicas=50

for dir in XXXsystemXXX ;
do
    echo $dir
    cd $dir 2> /dev/null

    version=$(echo $dir | sed "s/prodruns_//g")
    echo $version

    mkdir -p ${wdir}

    if [[ -e ${wdir}/output.txt ]];
    then
	echo "Analysis ongoing in ${dir}. Calculation aborted!"
	continue
    fi

    cd $wdir

    outdir=contact_maps_at_fixed_timestep
    mkdir -p ${outdir}

    gzfile=$(ls -1 *6_6667*gz)
    echo $gzfile

    rm -fr _tmp_fixed_*
    cat <( zcat ${gzfile} ) | awk '{if(NF==5 && $5<=20000000){f="_tmp_fixed_"$5;print $0 >> f}}'
    
    for rc in 6.6667 ; 
    do
	dcutoff=$(awk -v rc=$rc 'BEGIN{print int(rc*30)}')
	
	for res in 30kb ;
	do

	    outfile=${inputHiC}/cumulative_norm_HiC_matrices_4_dna_at_${res}.txt
	    
	    awk '{if($1=="#" && $2 == "MASKED"){for(i=3;i<=NF;i++) rem[$3]=$3; next}; if($1!="#"){n++; for(i=1;i<=n;i++){if(n in rem || i in rem){print n-1,i-1,"nan"}else{print n-1,i-1,$i}}}}' <( awk '{if(NF!=0) print $0}' ${outfile} ) > _exp
	    
	    # Produce genome-wide model contact map
	    resolution=$(echo ${res} | sed "s/kb/000/g")
	    
	    # Bins of the NOR to remove
	    MappingFile=bin_mapping_at_${res}_chr4.txt
	    
	    echo $rc ${resolution}
	    size=$(tail -1 ${scriptsdir}/Utilities/${MappingFile} | awk '{print $4+1}')
	    format=mat

	    #for timestep in $(seq 0 500000 1000000);
	    for timestep in $(seq 1000000 500000 20000000);
	    do
		outfile=${outdir}/contact_matrix_a_thaliana_chr_4_at_${rc}sigma_at_${res}_at_${timestep}_${version}.${format}
		
		if [[ ! -e ${outfile} ]];
		then
		    echo ${outfile}
		    
		    if [[ ${format} == "mat" ]];
		    then

			if [[ -e _tmp_fixed_${timestep} ]];
			then

			    awk -v rc=${rc} '{if(NF==4){b[$1-1]=$4;}; if(NF==5 && $4<=(rc*rc)){ind1=b[$2] ; ind2=b[$3] ; print ind1"_"ind2}}' ${scriptsdir}/Utilities/${MappingFile} _tmp_fixed_${timestep} | awk '{h[$1]++}END{for(i in h) print i,h[i]}' | sed "s/_/ /g" | awk -v s=${size} 'BEGIN{for(i=0;i<s;i++)for(j=0;j<s;j++) m[i,j]=0}{if(NF==3){i=int($1);j=int($2);m[i,j]=$3;m[j,i]=$3;}}END{for(i=0;i<s;i++){for(j=0;j<s;j++){printf("%s ",m[i,j])};printf("\n")}}' > ${outfile}

			    awk '{for(i=1;i<=NF;i++){print NR-1,i-1,$i}}' ${outfile} > _model
			    corr=$(python ${scriptsdir}/Utilities/compute_correlation_between_models_and_HiC_matrices_at_fixed_timestep.py $size $dcutoff $timestep 2> /dev/null)
			    mv plot.png ${outfile%.mat}.png
			    #head ${outfile}
			    bash ~/matrix_row_column_count.sh ${outfile}
			    
			fi
		    
			#tail -3 ${outfile} #| awk '{print $(NF-2),$(NF-1),$NF}'
			#rm -fr _tmp_fixed_${timestep}
		    fi
		fi
		outfile=${outdir}/contact_matrix_a_thaliana_chr_4_at_${rc}sigma_at_${res}_at_${timestep}_cumulative_${version}.${format}
		if [[ ! -e ${outfile} ]];
		then		
		    if [[ ${timestep} -eq 0 ]];
		    then			
			infile=${outdir}/contact_matrix_a_thaliana_chr_4_at_${rc}sigma_at_${res}_at_${timestep}_${version}.${format}
			cp ${infile} ${outfile}
		    fi
		    if [[ ${timestep} -gt 0 ]];
		    then			
			t=$((${timestep}-500000))
			infile1=${outdir}/contact_matrix_a_thaliana_chr_4_at_${rc}sigma_at_${res}_at_${timestep}_${version}.${format}
			infile2=${outdir}/contact_matrix_a_thaliana_chr_4_at_${rc}sigma_at_${res}_at_${t}_cumulative_${version}.${format}
			
			#awk -v s=${size} 'BEGIN{for(i=0;i<s;i++)for(j=0;j<s;j++)m[i,j]=0.0}{for(i=1;i<=s;i++){m[i-1,j-1]=$i+$(i+s)}}END{for(i=0;i<s;i++){for(j=0;j<s;j++){printf("%s ",m[i,j])};printf("\n")}}' <(paste ${infile1} ${infile2}) > ${outfile}
			awk -v s=${size} 'BEGIN{for(i=0;i<s;i++)for(j=0;j<s;j++)m[i,j]=0.0}{for(i=1;i<=s;i++){m[i-1,NR-1]=$i+$(i+s)}}END{for(i=0;i<s;i++){for(j=0;j<s;j++){printf("%s ",m[i,j])};printf("\n")}}' <(paste ${infile1} ${infile2}) > ${outfile}
		    fi
		    awk '{for(i=1;i<=NF;i++){print NR-1,i-1,$i}}' ${outfile} > _model
		    corr=$(python ${scriptsdir}/Utilities/compute_correlation_between_models_and_HiC_matrices_at_fixed_timestep.py $size $dcutoff $timestep 2> /dev/null)
		    mv plot.png ${outfile%.mat}.png
		    
		fi
	    done # Close cycle over ${timestep}
	    
	done # Close cycle over ${res}
    done # Close cycle over ${rc}
    rm -fr _tmp_fixed_* _model _exp
    
    cd .. # Exit ${wdir}
    cd .. # Exit ${dir}

done # Close cycle over ${dir}
echo "Analysis DONE!"
