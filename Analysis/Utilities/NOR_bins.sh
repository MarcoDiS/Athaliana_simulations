scriptsdir="."

for resolution in 3000 ;
do
    
    cat <(grep "A" ${scriptsdir}/INPUT_PARAMETERS_COMPLETE.txt | awk -v r=${resolution} '{print $1,int($3*(3000/r))+1}' | sed "s/A//g") <(grep A ${scriptsdir}/NORs_input.txt | awk -v r=${resolution} '{print $1,$2,int($3/r)+1,int($4/r)+1}' | sed "s/A//g") | awk 'BEGIN{offset=0}{if(NF==2){o[$1]=offset;offset+=$2}else{print $1,$2,$3+o[$2],$4+o[$2]}}' > ${scriptsdir}/NOR_at_${resolution}bp.txt

done
