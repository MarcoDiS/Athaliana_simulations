scriptsdir=/scratch/devel/mstefano/2018_05_04_Project_A_thaliana_physics_models/WT_simulations/Initial_conformations/generate_initial_conformation_v_chromosomes/scripts_and_programs/

# Obs-over-exp is compute on the full (genome-wide) matrix removing the low-coverage bins.
# Expected is the mean number of contact at that genomic distance over all the chromosomes.
# Write the matrices in a .tab format reporting the chromosome
for resolution in 30kb ; #150kb 30kb 15kb ;
do

    echo $resolution
    infile=cumulative_norm_HiC_matrices_full_at_${resolution}.mat
    infilecontacts=contacts_vs_gendist_a_thaliana_norm_full_${resolution}.txt
    outfile=cumulative_norm_HiC_matrices_full_at_${resolution}.tab

    ls -1 ${infile} ${infilecontacts} ${outfile}

    bash ~/matrix_row_column_count.sh ${infile}
    size=$(bash ~/matrix_row_column_count.sh ${infile} | head -1 | awk '{print $NF}')
    echo "Size $size"
    tail -1 ${scriptsdir}/bin_mapping_at_${resolution}.txt

    res=$(echo $resolution | sed "s/kb/000/g")
    grep "#" ${infile} | grep -v "MASKED" > ${outfile}
    grep "#" ${infile} | grep "MASKED" | awk '{printf("%s %s ",$1,$2);for(i=3;i<=NF;i++) printf("%d ",$i-1); printf("\n")}' >> ${outfile}
    grep "#" ${infile} | grep "MASKED" | awk '{for(i=3;i<=NF;i++) print $i-1}' > _masked

    # OK intra- or inter-chromosome contacts between well defined bins
    # NA inter-chromosome contacts between well defined bins
    # MASKED bad columns intra- or inter-chromosome
    echo "#chri parti chrj partj Obs Flag" >> ${outfile}
    cat _masked <(awk '{print $1,$2}' ${infilecontacts}) <(awk '{if(NF==4) print $4,$2,"N"}' ${scriptsdir}/bin_mapping_at_${resolution}.txt) <(awk '{for(i=1;i<=NF;i++){if(i<=NR)print NR-1,i-1,$i,"N"}}' <(grep -v "#" ${infile})) | awk -v r=${res} '{if(NF==1){m[$1]=$1}; if(NF==2){e[$1/r]=$2;}; if(NF==3){chr[$1]=$2}; if(NF==4){c1=chr[$1]; c2=chr[$2]; corrected="NA"; if(c1==c2){if($2<=$1){d=$1-$2}else{d=$2-$1}; if(d!=0){corrected="OK"}else{corrected="OK"}}; if($1 in m){corrected="MASKED"}; if($2 in m){corrected="MASKED"}; print c1,$1,c2,$2,$3,corrected}}' >> ${outfile}
    rm _masked
done
