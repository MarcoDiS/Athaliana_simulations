min_timestep=XXXmintimestepXXX

scriptsdir=../../Analysis/
input=${scriptsdir}/Utilities/contact_maps/

system=contact_map
indir=input_files

#maxrc=14.0
#maxrct=14_0
maxrc=6.6667
maxrct=6_6667

nreplicas=10

for dir in XXXsystemXXX ;
do
    echo $dir
    cd $dir #2> /dev/null
    

    version=$(echo $dir | sed "s/prodruns_//g")
    echo $version

    mkdir -p ${system}

    echo "Calculating the contacts within ${maxrc}sigma distance"
    for replica in XXXreplicaXXX ; 
    do

	# Getting the files for calculation
	replicadir=$(ls -1 | grep replica_${replica} | awk '{print $1}')
	cd $replicadir
	
	echo "Getting the name of the XYZ file"
	dir=$PWD
	check=$(ls -1 c* | grep XYZ | wc -l | awk '{print $1}' );
	if [[ $check -ne 1 ]] ;
	then
	    echo "ERROR: There are 2 XYZ files in replica_${replica}: ${check}"
	    cd .. # Exit replica
	    exit
	fi
	XYZ_file=$(ls -1 c* | grep XYZ )
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
        #chrs10_6_prodruns_replica_1
	name=$(echo ${XYZ_file%.XYZ})
	echo $name

	awk -v name=$name -v r=${replica} -v N=$N 'BEGIN{c=0;}{if(NR%N==1){next};if(NR%N==2){c=$1; n=name"_"c".XYZ"; print "ITEM: NUMBER OF ATOMS" > n; print $0 > n; next};print $0 > n;if(NR%N==0){close(n)}}' $XYZ_file 

	if [[ ${min_timestep} -gt 0 ]];
	then
	    echo "Min timestep: ${min_timestep}"
	    for timestep in $(seq 0 ${delta_timestep} ${min_timestep} );
	    do
		if [[ -e ${name}_${timestep}.XYZ ]];
		then
		    rm -fvr ${name}_${timestep}.XYZ
		fi
	    done
	fi
	for timestep in $(seq 20500000 ${delta_timestep} 40000000 );
	do
	    if [[ -e ${name}_${timestep}.XYZ ]];
	    then
		rm ${name}_${timestep}.XYZ
	    fi
	done	

	mv -v ./${name}_*.XYZ ..
	cd ..

	mkdir -p $system
	echo $system
	
	cd $system
	rm -fvr ./*_replica_${replica}_*XYZ
	mv ../*_replica_${replica}_*XYZ .    

        # Check whether the calculation is done
	check1=9 #$(ls -1 *gz 2> /dev/null | wc -l)
	gzfile=distances_model_a_thaliana_replica_${replica}_below_${maxrct}sigma_at_3kb_txt.tar.gz
	if [[ ! -e ${gzfile} ]];
	then
	    contactfile=distances_model_a_thaliana_replica_${replica}_below_${maxrc}sigma_at_3kb.txt	    
	    # If it is not done, let's try to start it!
	    if [[ -e ${outfile} ]];
	    then
		echo "The calculation is already ongoing. Aborting!"
		exit
	    fi
	    tmpdir=_tmp_replica_${replica}
	    if [[ -e ${tmpdir}/output.txt ]];
	    then
		echo "Analysis ongoing in ${PWD}/${tmpdir}. Calculation aborted!"
		exit
	    fi
	    # Since it is not started, let's do it!
	    mkdir -p ${indir}	
	    mkdir -p ${tmpdir}
	    cd ${tmpdir}
	    mv ../*_replica_${replica}_*XYZ .

	    echo "#Replica i j cij time" > ${contactfile}
	    echo "replica_${replica}"

	    for file in $(ls -1 *_replica_${replica}_*XYZ);
	    do

		ls -1 ${file} > DATA_FILES_INPUT.txt	
		nparticles=$(tail -1 ${file} | awk '{print $1}')
		awk -v np=${nparticles} '{if($3<=np) print $0}' $input/mapping_particles_single_chromosome.txt > mapping_particles.txt
		    
		time=$(echo ${file} | sed -e "s/\.XYZ/ /g" -e "s/_/ /g" | awk '{print $NF}')
		echo $(cat D*txt) $nparticles $time
		    
                # Resolution 3kb -> 1 particles we could do the coarsening later
		#echo $outfile
		${scriptsdir}/contact_maps/contacts_a_thaliana -p ${nparticles} -r 1 -d ${maxrc} -k 2
		cat output.txt | awk -v r=${replica} -v t=${time} '{print r,$0,t}' | grep -v XYZ >> ../${contactfile} 
		rm -fr output.???
		mv ${file} ../${indir}
	    done # Close cycle over ${file}
	    cd .. # Exit ${tmpdir}
	    bash ~/create_tarball.sh ${contactfile}
	    rm -fr ${contactfile} ${tmpdir}
	fi # Close condition over ${gzfile}
	cd .. # Exit $system
    done # Close cycle over ${replica}

    cd $system

    check2=$(ls -1 *gz 2> /dev/null | wc -l)
    if [[ $check1 -eq $((${nreplicas}-1)) && $check2 -eq $((${nreplicas})) ]];
    then
	for rc in 6.6667 ;
	do

	    echo "dcutoff = ${rc}"
	    cat <( zcat distances_model_a_thaliana_replica_*_below_${maxrct}sigma_at_3kb_txt.tar.gz ) | awk -v rc=${rc} -v s=${size} '{if($4<=(rc*rc) && $5>=1000000 && $5<=20000000){print $0}}' > _tmp

	    for res in 30kb ; #150kb 30kb 15kb ;
	    do

		# Produce genome-wide model contact map
		resolution=$(echo ${res} | sed "s/kb/000/g")

		# Bins of the NOR to remove
		MappingFile=bin_mapping_at_${res}_chr4.txt
		
		echo $rc ${resolution}
		if [[ ! -e ${scriptsdir}/Utilities/${MappingFile} ]];
		then
		    echo "ERROR! File ${scriptsdir}/Utilities/${MappingFile} missing"
		    echo "Create it first using ${scriptsdir}/Utilities/bin_mapping.sh"
		    continue
		fi
		size=$(tail -1 ${scriptsdir}/Utilities/${MappingFile} | awk '{print $4+1}')
		#format=mat
		format=tab
		outfile=contact_matrix_a_thaliana_at_${rc}sigma_at_${res}_${version}.${format}
		outfiledecay=contacts_vs_gendist_a_thaliana_at_${rc}sigma_at_${res}_${version}.txt
		outfilenormdecay=norm_contacts_vs_gendist_a_thaliana_at_${rc}sigma_at_${res}_${version}.txt
		echo $size

		if [[ ! -e ${outfile} ]];
		then
		    echo ${outfile}
		    
		    if [[ ${format} == "mat" ]];
		    then
			# Mapping particles to bins not considering the NORs  
			cat <( zcat distances_model_a_thaliana_below_${maxrct}sigma_at_3kb_txt.tar.gz ) | awk -v rc=${rc} '{if($4<=(rc*rc)){print $0}}' > _tmp

			cat ${scriptsdir}/Utilities/${MappingFile} _tmp | awk '{if(NF==4){b[$1-1]=$4;}; if(NF==5){ind1=b[$2] ; ind2=b[$3] ; print ind1"_"ind2}}' | awk '{h[$1]++}END{for(i in h) print i,h[i]}' | sed "s/_/ /g" | awk -v s=${size} 'BEGIN{for(i=0;i<s;i++)for(j=0;j<s;j++) m[i,j]=0}{if(NF==3){i=int($1);j=int($2);m[i,j]=$3;m[j,i]=$3}}END{for(i=0;i<s;i++){for(j=0;j<s;j++){printf("%s ",m[i,j])};printf("\n")}}' > ${outfile}

			bash ~/matrix_row_column_count.sh ${outfile}

			rm -fr _tmp


			if [[ ! -e ${outfile%.mat}.pdf ]];
			then
			    cp ${outfile} raw_matrix.txt
			    sz=$(wc -l raw_matrix.txt | awk '{print $1}')
			    echo $sz
			    python ${scriptsdir}/Utilities/plot_contact_matrix.py ${sz}
			    mv raw_matrix.pdf ${outfile%.mat}.pdf
			    ls -lrth ${outfile%.mat}.pdf
			    rm -fr raw_matrix.txt
			fi		

		    fi

		    if [[ ${format} == "tab" ]];
		    then

			awk -v s=${size} 'BEGIN{for(i=0;i<s;i++)for(j=0;j<=i;j++)print 11,i,j,0,0}' > _tmp_add

			#echo "# CRM 1_dna 30427671" >  ${outfile}
			#echo "# CRM 2_dna 19698289" >> ${outfile}
			#echo "# CRM 3_dna 23459830" >> ${outfile}
			echo "# CRM 4_dna 18585056" >> ${outfile}
			#echo "# CRM 5_dna 26975502" >> ${outfile}
			echo "#chri parti chrj partj Obs Obs/Exp" >> ${outfile}

			cat <(cat ${scriptsdir}/${MappingFile} _tmp | awk '{if(NF==4){b[$1-1]=$4; c[$1-1]=$2}; if(NF==5){if($2 in b) if($3 in b){ind1=b[$2] ; chr1=c[$2] ; ind2=b[$3] ; chr2=c[$3] ; if(ind2<=ind1){print chr1"_"ind1"_"chr2"_"ind2}else{print chr2"_"ind2"_"chr1"_"ind1}}}}') <(cat ${scriptsdir}/${MappingFile} _tmp_add | awk '{if(NF==4){b[$4]=$4; c[$4]=$2}; if(NF==5){if($2 in b) if($3 in b){ind1=b[$2] ; chr1=c[$2] ; ind2=b[$3] ; chr2=c[$3] ; if(ind2<=ind1){print chr1"_"ind1"_"chr2"_"ind2}else{print chr2"_"ind2"_"chr1"_"ind1}}}}') | awk '{h[$1]++}END{for(i in h) print i,h[i]-1}' | sed "s/_/ /g" | sort -k 2n,2n -k 4n,4n > _tmp_matrix

			awk '{if($1==$3){d=$2-$4;print d,$5}}' _tmp_matrix | awk '{e[$1]+=$2;e2[$1]+=($2*$2);c[$1]++}END{for(i in e){avg=e[i]/c[i];avg2=e2[i]/c[i];stddev=sqrt(avg2-avg*avg);print i,avg,stddev}}' > _tmp_decay
			cat _tmp_decay _tmp_matrix | awk '{if(NF==3){w[$1]=$2};if(NF==5){corrected="NA";if($1==$3){d=$2-$4; if(w[d]!=0){corrected=$5/w[d]}else{corrected=$5}}; print $0,corrected}}' >> ${outfile}

			awk -v r=${resolution} '{c[$1*r]=$2;e[$1*r]=$3;f=1}END{for(i in e) print i,c[i]/f,e[i]/f}' _tmp_decay | sort -k 1n > ${outfiledecay}
			awk -v r=${resolution} '{c[$1*r]=$2;e[$1*r]=$3;if($1*r==300000)f=$2}END{for(i in e) print i,c[i]/f,e[i]/f}' _tmp_decay | sort -k 1n > ${outfilenormdecay}

			rm -fvr _tmp_matrix _tmp_decay _tmp_add
		    fi
		    
		    tail -3 ${outfile} #| awk '{print $(NF-2),$(NF-1),$NF}'
		    
		fi

	    done
	    rm -fvr _tmp
	done
    fi

    cd .. # Exit ${system}
    cd .. # Exit ${dir}

done
echo "Analysis DONE!"
