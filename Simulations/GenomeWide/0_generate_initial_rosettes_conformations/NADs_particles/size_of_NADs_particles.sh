R=46.6666
nNADs_particles=$(cat NADs_particles_NORs_only.txt | sed "s/*/ /g" | awk '{s+=$4-$3+1}END{print s}')
echo $R $nNADs_particles

awk -v nN=${nNADs_particles} -v r=${R} 'BEGIN{printf("Target radius for NADs particles (rescaled for optimal FCC packing) = %f factor %f\n",(r*r*r*0.74/nN)^(1/3), 3.14159265359/(3.*sqrt(2)))}'
awk -v nN=${nNADs_particles} -v r=${R} 'BEGIN{printf("Target radius for NADs particles (rescaled for random  RCP packing) = %f factor %f\n",(r*r*r*0.64/nN)^(1/3), 0.64)}'
awk -v nN=${nNADs_particles} -v r=${R} 'BEGIN{printf("Target radius for NADs particles                                    = %.1f\n",(r*r*r/nN)^(1/3))}'