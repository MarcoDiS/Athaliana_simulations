
scriptsdir=../../Analysis/
n_particles=84502
s=$1 #State

awk '{print $1,$2}' ${scriptsdir}/Utilities/chromatin_states_genome_wide.txt | grep -v "#" | sort -k 2n | uniq > _list_of_particles_${s}
pi=3.14159265359

#3  - Associate the particle position to the corresponding x,rho coordinates. The grid is a square of size 2500.0x2500.0nm and the resolution is 60nm (2.0). Creating a mask to filter the particles
awk -v n_particles=$n_particles 'BEGIN{mask[n_particles]; tag[n_particles];for(i = 0;i < n_particles;i++){mask[i]=0;tag[i]="NULL"}}{mask[$1-1]++;tag[$1-1]=$2}END{for(i = 0;i < n_particles;i++) print mask[i],tag[i],i+1}' _list_of_particles_${s} > _mask_${s}
wc -l _mask_${s}

i=1
wc -l DATA_FILES_INPUT.DAT
nfiles=$(wc -l DATA_FILES_INPUT.DAT | awk '{print $1}')
voxelfile=radial_position_all_file.txt
if [[ $s == "All" ]];
then 
    echo "I'm here!"
    if [[ ! -e ${voxelfile} ]];
    then
	echo "#RadPos state replica timestep" > ${voxelfile}
	for file in $(cat DATA_FILES_INPUT.DAT) ;    
	do
	    echo "Analysing file "$i" over "$nfiles" files"
	    i=$(($i+1))
	    echo $file
	    #chrs10_6_prodruns_replica_50_16000000
	    replica=$(echo $file | sed "s/_/ /g" | awk '{print $5}')
	    time=$(echo $file | sed -e "s/\.XYZ//g" -e "s/_/ /g" | awk '{print $6}')
	    xcom=$(grep sphere ../../replica_${replica}/_tmp.lmp | awk '{print $4}')
	    ycom=$(grep sphere ../../replica_${replica}/_tmp.lmp | awk '{print $5}')
	    zcom=$(grep sphere ../../replica_${replica}/_tmp.lmp | awk '{print $6}')
	    echo ${replica} ${xcom} ${ycom} ${zcom} 
	    awk -v x=${xcom} -v y=${ycom} -v z=${zcom} '{if(NR>9){print $1,$3-x,$4-y,$5-z}}' $file  > input.dat
	
	    #Creating the list of coordinates
	    awk -v r=${replica} -v t=${time} '{print $1,r,t,sqrt($2*$2+$3*$3+$4*$4)}' input.dat > _coord_${s}
	    
	    #Filtering the coordintes applying the mask. Output: particle index, discretized x and y
	    paste _mask_${s} _coord_${s} | awk '{if($3==$4) print $7,$2,$5,$6}' >> ${voxelfile}
	done
	rm -fvr input.dat
    fi
    rm -fvr _*${s}
    echo "DONE general calculation!"
fi

if [[ $s == "All" ]];
then
    exit
fi

for state in $s ;
do
    echo $state
    
    if [[ ${state} == "Active" ]];
    then
	color=red
    fi
    if [[ ${state} == "Heterochromatin" ]];
    then
	color=green
    fi
    if [[ ${state} == "Polycomb_like" ]];
    then
	color=blue
    fi
    if [[ ${state} == "NORs" ]];
    then
	color=gray
    fi
    if [[ ${state} == "Telomere" ]];
    then
	color=gray
    fi
    if [[ ${state} == "Undetermined" ]];
    then
	color=black
    fi

    # Equal width bins
    for Binwidth in 125 250 500 ; 
    do
	
	outfile=${state}_histogram_at_${Binwidth}nm_radial_position.txt
	if [[ ! -s ${outfile} ]];
	then
	    
	    if [[ ! -s radial_position_${state}.txt ]];
	    then
		grep -w ${state} ${voxelfile} > radial_position_${state}.txt
	    fi
	    nreplicas=$(awk '{print $(NF-1)}' radial_position_${state}.txt | uniq | sort -k 1n | tail -1 | awk '{print $1}')    
	    echo $nreplicas    

	    awk -v pi=$pi -v nr=${nreplicas} -v bw=${Binwidth} '{if(NF==4){bin=int(($1*30)/bw); h[$(NF-1),bin]++; cnt[$(NF-1)]++;}; if(NF==3){bin=int(($1*30)/bw); h[$NF,bin]++; cnt[$NF]++;}}END{for(i=1;i<=nr;i++){for(j=0;j<=((2500-bw)/bw);j++){R=(j+1)*bw; r=(j+0)*bw; V=4./3.*pi*(R*R*R-r*r*r); if(cnt[i]==0){continue}else{print j*bw,h[i,j]/(cnt[i]*V),i}}}}' radial_position_${state}.txt | sort -k 1n > ${outfile}

	    rm -fr _a_${state}
	    for bin in $(awk '{print $1}' ${outfile} | sort | uniq);
	    do
		#echo $bin
		if [[ ${nreplicas} -ne 1 ]];
		then
		    awk -v b=${bin} '{if($1==b) print $2}' ${outfile} | ~/howmuch.csh | awk -v b=${bin} '{print b,$10,$14}' >> _a_${state}
		fi
		if [[ ${nreplicas} -eq 1 ]];
		then
		    awk -v b=${bin} '{if($1==b) print $0}' ${outfile} >> _a_${state}
		fi
	    done
	    
	    sort -k 1n _a_${state} > ${outfile}
	    cat ${outfile}	    
	fi
    done
done
rm -fvr _*${s}
