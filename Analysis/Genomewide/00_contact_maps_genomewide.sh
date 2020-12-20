min_timestep=XXXmintimestepXXX

scriptsdir=../../Analysis/
input=${scriptsdir}/Utilities/contact_maps/

wdir=contact_map
indir=input_files
format=tab

nreplicas=50

for dir in XXXsystemXXX ;
do
    echo $dir
    cd $dir 2> /dev/null

    version=$(echo $dir | sed "s/prodruns_//g")
    echo $version


    mkdir -p ${wdir}

    if [[ -e ${wdir}/outoput.txt ]];
    then
	echo "Analysis ongoing in ${dir}. Calculation aborted!"
	continue
    fi

    rc=XXXrcXXX
    rcstr=$(echo ${rc} | sed "s,\.,_,g")
    echo "Calculating the distance matrix"

    for replica in XXXreplicaXXX ;
    do
	echo "Compute contacts per replica ${replica}"

	gzfile=distances_model_a_thaliana_replica_${replica}_below_${rcstr}sigma_at_3kb_txt.tar.gz
	outfile=distances_model_a_thaliana_replica_${replica}_below_${rc}sigma_at_3kb.txt

	if [[ ! -e ${wdir}/${gzfile} ]];
	then
	    if [[ -e ${wdir}/distances_model_a_thaliana_replica_${replica}_below_${rcstr}sigma_at_3kb.txt ]];
	    then
		echo "The calculation is already ongoing. Aborting!"
		exit
	    fi	    
	    replicadir=$(ls -1 | grep replica_${replica} | awk '{print $1}')
	    cd $replicadir
	    
	    echo "Getting the name of the XYZ file"
	    dir=$PWD
	    check=$(ls -1 c* | grep XYZ | grep -v part | wc -l | awk '{print $1}' );
	    if [[ $check -ne 1 ]] ;
	    then
		echo "There are 2 XYZ files in replica_${replica}: ${check}"
		cd .. # Exit replica
		continue
	    fi
	    XYZ_file=$(ls -1 c* | grep XYZ |  grep -v part )
	    echo $XYZ_file
	    
	    echo "Getting the number of particles"
	    L=$(head -4 $XYZ_file | tail -1 | awk '{print $1}')
	    N=$(($L+9))
	    echo "Number of particles "$L $N
	    
	    echo "Storing the last timestep and the delta between timesteps"
	    last_timestep=$(tail -${N} $XYZ_file | head -2 | tail -1 | awk '{print $1}')
	    delta_timestep=$(tail -$(($N*2)) $XYZ_file | awk -v N=$N '{if(NR==2){f=$1;next};if(NR%N==2){print $1-f}}')
	    echo "Last and delta timestep "$last_timestep $delta_timestep
	   
	    echo "Splitting the XYZ file"
	    name=$(echo ${XYZ_file%.XYZ})
	    echo $name
	    
	    awk -v name=$name -v r=${replica} -v N=$N 'BEGIN{c=0;}{if(NR%N==1){next};if(NR%N==2){c=$1; n=name"_"c".XYZ"; print "ITEM: NUMBER OF ATOMS" > n; print $0 > n; next};print $0 > n;if(NR%N==0){close(n)}}' $XYZ_file 
	    
	    for timestep in $(seq 0 ${delta_timestep} ${min_timestep} );
	    do
		if [[ -e ${name}_${timestep}.XYZ ]];
		then
		    rm ${name}_${timestep}.XYZ
		fi
	    done
	    
	    mv -v ./${name}_*0000*XYZ ..
	    cd ..

	    echo $wdir # Change directory from replica_${replica} to ${wdir}
	    mkdir -p $wdir
	    
	    cd $wdir
	    rm -fvr ./*_replica_${replica}_*XYZ
	    mv ../*_replica_${replica}_*XYZ .

	    mkdir -p ${indir}	

	    for rc in ${rc} ;
	    do
		# Enter in a tmp directory to avoid files overwriting when computing the contact
		# for different replicas at the same time
		tmpdir=_tmp_replica_${replica}
		mkdir -p ${tmpdir}
		cd ${tmpdir}
		mv ../*_replica_${replica}_*XYZ .

		echo "#Replica i j cij time" > ../${outfile}
		
		echo "replica_${replica}"		
		for file in $(ls -1 *_replica_${replica}_*XYZ);
		do
		    
		    ls -1 ${file} > DATA_FILES_INPUT.txt	
		    nparticles=$(tail -1 ${file} | awk '{print $1}')
		    awk -v np=${nparticles} '{if($3<=np) print $0}' $input/mapping_particles_entire_genome.txt > mapping_particles.txt
		    
		    time=$(echo ${file} | sed -e "s/\.XYZ/ /g" -e "s/_/ /g" | awk '{print $NF}')
		    echo $(cat D*txt) $nparticles $time
		    
                    # Resolution 3kb -> 1 particles we could do the coarsening later
		    ${scriptsdir}/contact_maps/contacts_a_thaliana -p ${nparticles} -r 1 -d ${rc} -k 2
		    cat output.txt | awk -v r=${replica} -v t=${time} '{print r,$0,t}' | grep -v XYZ >> ../${outfile} 
		    rm -fr output.???
		    mv ${file} ../${indir}
		done # Close cycle over ${file}
		cd .. # Exit from ${tmpdir}
		rm -fr ${tmpdir}
		bash ~/create_tarball.sh ${outfile}
		rm -fr ${outfile}
	    done # Close cycle over ${replica}
	    cd .. # Exit $wdir
	fi

	echo "Computing contact maps per replica ${replica}"
	cd $wdir
	
	if [[ -e ${gzfile} ]];
	then	
	    for rc in ${rc} ;
	    do

		# Check whether there are already all the contact maps
		check=0
		for res in 30kb ;
		do
		    outfile=contact_matrix_a_thaliana_at_${rc}sigma_at_${res}_replica_${replica}_${version}.${format}
		    if [[ ! -e ${outfile} ]];
		    then
			check=1
		    fi
		done
		if [[ $check -eq 0 ]];
		then
		    continue
		fi
		
		echo "dcutoff ${rc}"
		tmpadd=_tmp_add_${replica}
		tmpmatrix=_tmp_matrix_${replica}
		tmpdecay=_tmp_decay_${replica}
		tmpreplica=_tmp_replica_${replica}
		cat <( zcat ${gzfile} ) | awk -v rc=${rc} -v s=${size} '{if($4<=(rc*rc)){print $0}}' > ${tmpreplica}
		
		for res in 30kb ; #150kb #15kb ;
		do
		    
		    # Produce genome-wide model contact map
		    resolution=$(echo ${res} | sed "s/kb/000/g")
		    
		    # Bins of the NOR to remove
		    NORfile=NORs_output_at_${res}.txt
		    
		    echo $rc ${resolution}
		    if [[ ! -e ${scriptsdir}/bin_mapping_at_${res}.txt ]];
		    then
			echo "ERROR! File ${scriptsdir}/bin_mapping_at_${res}.txt missing"
			echo "Create it first using ${scriptsdir}/bin_mapping.sh"
			continue
		    fi
		    size=$(tail -1 ${scriptsdir}/bin_mapping_at_${res}.txt | awk '{print $4+1}')

		    outfile=contact_matrix_a_thaliana_at_${rc}sigma_at_${res}_replica_${replica}_${version}.${format}
		    outfiledecay=contacts_vs_gendist_a_thaliana_at_${rc}sigma_at_${res}_replica_${replica}_${version}.txt
		    outfilenormdecay=norm_contacts_vs_gendist_a_thaliana_at_${rc}sigma_at_${res}_replica_${replica}_${version}.txt
		    echo $size
		    
		    if [[ ! -e ${outfile} ]];
		    then
			echo ${outfile}
			
			if [[ ${format} == "tab" ]];
			then

			    awk -v s=${size} 'BEGIN{for(i=0;i<s;i++)for(j=0;j<=i;j++)print 51,i,j,0,0}' > ${tmpadd}
			    
			    echo "# CRM 1_dna 30427671" >  ${outfile}
			    echo "# CRM 2_dna 19698289" >> ${outfile}
			    echo "# CRM 3_dna 23459830" >> ${outfile}
			    echo "# CRM 4_dna 18585056" >> ${outfile}
			    echo "# CRM 5_dna 26975502" >> ${outfile}
			    echo "#chri parti chrj partj Obs Obs/Exp" >> ${outfile}
			    
			    cat <(cat ${scriptsdir}/bin_mapping_at_${res}.txt ${tmpreplica} | awk '{if(NF==4){b[$1-1]=$4; c[$1-1]=$2}; if(NF==5){if($2 in b) if($3 in b){ind1=b[$2] ; chr1=c[$2] ; ind2=b[$3] ; chr2=c[$3] ; if(ind2<=ind1){print chr1"_"ind1"_"chr2"_"ind2}else{print chr2"_"ind2"_"chr1"_"ind1}}}}') <(cat ${scriptsdir}/bin_mapping_at_${res}.txt ${tmpadd} | awk '{if(NF==4){b[$4]=$4; c[$4]=$2}; if(NF==5){if($2 in b) if($3 in b){ind1=b[$2] ; chr1=c[$2] ; ind2=b[$3] ; chr2=c[$3] ; if(ind2<=ind1){print chr1"_"ind1"_"chr2"_"ind2}else{print chr2"_"ind2"_"chr1"_"ind1}}}}') | awk '{h[$1]++}END{for(i in h) print i,h[i]-1}' | sed "s/_/ /g" | sort -k 2n,2n -k 4n,4n > ${tmpmatrix}

			    awk '{if($1==$3){d=$2-$4;print d,$5}}' ${tmpmatrix} | awk '{e[$1]+=$2;e2[$1]+=($2*$2);c[$1]++}END{for(i in e){avg=e[i]/c[i];avg2=e2[i]/c[i];stddev=sqrt(avg2-avg*avg);print i,avg,stddev}}' > ${tmpdecay}
			    cat ${tmpdecay} ${tmpmatrix} | awk '{if(NF==3){w[$1]=$2};if(NF==5){corrected="NA";if($1==$3){d=$2-$4; if(w[d]!=0){corrected=$5/w[d]}else{corrected=$5}}; print $0,corrected}}' >> ${outfile}
			    
			    awk -v r=${resolution} '{c[$1*r]=$2;e[$1*r]=$3;f=1}END{for(i in e) print i,c[i]/f,e[i]/f}' ${tmpdecay} | sort -k 1n > ${outfiledecay}
			    awk -v r=${resolution} '{c[$1*r]=$2;e[$1*r]=$3;if($1*r==300000)f=$2}END{for(i in e) print i,c[i]/f,e[i]/f}' ${tmpdecay} | sort -k 1n > ${outfilenormdecay}
			    
			    rm -fvr ${tmpmatrix} ${tmpdecay} ${tmpadd}
			fi
			
			tail -3 ${outfile} #| awk '{print $(NF-2),$(NF-1),$NF}'
			
		    fi
		    
		done
		rm -fvr ${tmpreplica}
	    done
	fi    
	cd .. # Exit ${wdir}
    done # Close cycle over {replica}

    if [[ ${replica} -ne 50 ]];
    then
	exit
    fi

    echo "Computing genome-wide contact maps"
    cd $wdir
	
    for rc in ${rc} ;
    do
	echo "dcutoff ${rc}"
	
	for res in 30kb ;
	do
	    resolution=$(echo ${res} | sed "s/kb/000/g")	    

	    # Produce genome-wide model contact map	    
	    outfile=contact_matrix_a_thaliana_at_${rc}sigma_at_${res}_${version}.${format}
	    outfiledecay=contacts_vs_gendist_a_thaliana_at_${rc}sigma_at_${res}_${version}.txt
	    outfilenormdecay=norm_contacts_vs_gendist_a_thaliana_at_${rc}sigma_at_${res}_${version}.txt
	    echo $size

	    if [[ -e ${outfile} ]];
	    then
		continue
	    fi

	    if [[ ! -e ${outfile} ]];
	    then
		echo ${outfile}
		
		echo "# CRM 1_dna 30427671" >  ${outfile}
		echo "# CRM 2_dna 19698289" >> ${outfile}
		echo "# CRM 3_dna 23459830" >> ${outfile}
		echo "# CRM 4_dna 18585056" >> ${outfile}
		echo "# CRM 5_dna 26975502" >> ${outfile}
		echo "#chri parti chrj partj Obs Obs/Exp" >> ${outfile}
		awk '{if(substr($1,1,1)=="#"){next}else{ind=$1"_"$2"_"$3"_"$4; h[ind]+=$5; cnt[ind]++}}END{for(i in h) print i,h[i],cnt[i]}' contact_matrix_a_thaliana_at_${rc}sigma_at_${res}_replica_*_${version}.${format} | sed "s/_/ /g" | sort -k 2n,2n -k 4n,4n > _tmp_matrix

		nrep=$(grep -v "#" _tmp_matrix | awk '{print $6}' | sort | uniq)
		if [[ ${nrep} -ne ${nreplicas} ]];
		then
		    echo "Number of replicas ${nrep} inconsistent with what expected ${nreplicas}"
		    exit
		fi

		awk '{if($1==$3){d=$2-$4;print d,$5}}' _tmp_matrix | awk '{e[$1]+=$2;e2[$1]+=($2*$2);c[$1]++}END{for(i in e){avg=e[i]/c[i];avg2=e2[i]/c[i];stddev=sqrt(avg2-avg*avg);print i,avg,stddev}}' > _tmp_decay
		cat _tmp_decay _tmp_matrix | awk '{if(NF==3){w[$1]=$2};if(NF==6){corrected="NA";if($1==$3){d=$2-$4; if(w[d]!=0){corrected=$5/w[d]}else{corrected=$5}}; print $1,$2,$3,$4,$5,corrected}}' >> ${outfile}
		
		awk -v r=${resolution} '{c[$1*r]=$2;e[$1*r]=$3;f=1}END{for(i in e) print i,c[i]/f,e[i]/f}' _tmp_decay | sort -k 1n > ${outfiledecay}
		awk -v r=${resolution} '{c[$1*r]=$2;e[$1*r]=$3;if($1*r==300000)f=$2}END{for(i in e) print i,c[i]/f,e[i]/f}' _tmp_decay | sort -k 1n > ${outfilenormdecay}
		
		rm -fvr _tmp_matrix _tmp_decay _tmp_add
	    fi
	    
	    tail -3 ${outfile} #| awk '{print $(NF-2),$(NF-1),$NF}'
	done
	rm -fvr _tmp
    done
    
    cd .. # Exit ${wdir}
    cd .. # Exit ${dir}
    
done
echo "Analysis DONE!"
