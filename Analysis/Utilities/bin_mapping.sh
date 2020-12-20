scriptsdir=.

for res in 3kb ;
do
    
    resolution=$(echo $res | sed "s/kb/000/g")

    cat ${scriptsdir}/INPUT_PARAMETERS_COMPLETE.txt | grep "A" | awk -v r=${resolution} 'BEGIN{p_off=0;b_off=0}{for(j=0;j<=1;j++){for(i=0;i<=($3-1);i++){print i+p_off+1,int(i*3000/r)+b_off};p_off=$3+p_off};b_off=int($3*3000/r)+b_off}' > ${scriptsdir}/bin_mapping_at_${res}.txt

    cat ${scriptsdir}/INPUT_PARAMETERS_COMPLETE.txt | grep "A" | awk -v r=${resolution} 'BEGIN{p_off=0;b_off=0}{for(j=0;j<=1;j++){for(i=1;i<=$3;i++){print $1,int(i*3000/r)+b_off};p_off=$3+p_off};b_off=int($3*3000/r)+b_off}' | sed "s/A//g" | sort -k 2n | uniq > ${scriptsdir}/chrs_mapping_at_${res}.txt

    awk 'BEGIN{n=0}{if(NF==4){s1[NR]=$3;s2[NR]=$4}else{for(i in s1) if(s1[i]<=$2 && $2<=s2[i]){print $1,$2; next} ; print $1,$2,n; n++}}' NOR_at_${res}.txt chrs_mapping_at_${res}.txt > _a ; mv _a chrs_mapping_at_${res}.txt
    exit

    awk '{if(NF>=3){if($3!="N")b[$2]=$3;c[$2]=$1}else{print $1,c[$2],$2,b[$2]}}' <( awk '{print $0,"N"}' chrs_mapping_at_${res}.txt) bin_mapping_at_${res}.txt > _a ; mv _a bin_mapping_at_${res}.txt

done
