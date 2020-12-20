cat <(cat ../0_generate_initial_rosettes_conformations/NADs_particles/NADs_particles_NORs_only.txt | sed "s/*/ /g" | awk '{print $3,$4}') all_mitotic_bonds.txt | awk '{if(NF==2){s1[NR]=$1;s2[NR]=$2}else{f=0;for(i in s1){if((s1[i]<=$7&&$7<=s2[i])||(s1[i]<=$6&&$6<=s2[i])){f=1}};if(f==0){print $0}}}' > noNADs_mitotic_bonds.txt

cat noNADs_mitotic_bonds.txt | awk '{print "un"$1,$2}' > unfix_noNADs_mitotic_bonds.txt

cat all_mitotic_bonds.txt | awk '{print "un"$1,$2}' > unfix_all_mitotic_bonds.txt

tail all_mitotic_bonds.txt unfix_all_mitotic_bonds.txt noNADs_mitotic_bonds.txt unfix_noNADs_mitotic_bonds.txt ; wc -l all_mitotic_bonds.txt unfix_all_mitotic_bonds.txt noNADs_mitotic_bonds.txt unfix_noNADs_mitotic_bonds.txt
