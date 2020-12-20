infile=NADs_particles_NORs_only.txt
outfile=normal_bonds.txt
outfileNAD=NADs_bonds.txt
outfileinter=interface_bonds.txt

cat <( cat ${infile} | sed "s/*/ /g" | awk '{print $3,$4}') bonds.txt | awk 'BEGIN{c=0}{if(NF==2){s1[c]=$1;s2[c]=$2;c++}else{c=0;for(i in s1){if((s1[i]<=$3&&$3<=s2[i])||(s1[i]<=$4&&$4<=s2[i])){c=1}};if(c==0){print $1,1,$3,$4}}}' > ${outfile}
cat <( cat ${infile} | sed "s/*/ /g" |awk '{print $3,$4}') bonds.txt | awk 'BEGIN{c=0}{if(NF==2){s1[c]=$1;s2[c]=$2;c++}else{for(i in s1){if((s1[i]<=$3&&$3<=s2[i])&&(s1[i]<=$4&&$4<=s2[i])){print $1,2,$3,$4}}}}' > ${outfileNAD}
cat <( cat ${infile} | sed "s/*/ /g" |awk '{print $3,$4}') bonds.txt | awk 'BEGIN{c=0}{if(NF==2){s1[c]=$1;s2[c]=$2;c++}else{for(i in s1){if((s1[i]<=$3&&$3<=s2[i])||(s1[i]<=$4&&$4<=s2[i])){print $1,2,$3,$4}}}}' > _tmp
diff _tmp ${outfileNAD} | awk '{if(NF==5)print $2,3,$4,$5}' > ${outfileinter}
echo "Inteface bonds "$(cat ${outfileinter})
rm -fr _tmp

cat <(echo) <(echo "Bonds") <(echo) > bonds_with_NADs_particles.txt
cat ${outfile} ${outfileNAD} ${outfileinter} | sort -k 1n | uniq >> bonds_with_NADs_particles.txt