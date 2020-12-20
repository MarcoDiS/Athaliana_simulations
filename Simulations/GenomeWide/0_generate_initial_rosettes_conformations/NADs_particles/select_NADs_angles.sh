infile=NADs_particles_NORs_only.txt
outfile=normal_angles.txt
outfileNAD=angles_with_NADs_particles.txt

cat <( cat ${infile} | sed "s/*/ /g" | awk '{print $3,$4}') angles.txt | awk 'BEGIN{c=0}{if(NF==2){s1[c]=$1;s2[c]=$2;c++}else{a=0;for(i in s1){if((s1[i]<=$3&&$3<=s2[i])||(s1[i]<=$4&&$4<=s2[i])||(s1[i]<=$5&&$5<=s2[i])){a=1}};if(a==0){print $1,1,$3,$4,$5}}}' | awk '{print NR,$2,$3,$4,$5}' > ${outfile}
wc -l ${outfile}

cat <(echo) <(echo "Angles") <(echo) > ${outfileNAD}
cat ${outfile} | sort -k 1n | awk '{print NR,$2,$3,$4,$5}' | uniq >> ${outfileNAD}