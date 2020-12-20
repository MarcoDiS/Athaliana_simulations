#!/bin/bash

maindir=${PWD}
scriptsdir=${PWD}/../scripts_and_programs/
genomicregionsfile=${PWD}/../scripts_and_programs/interesting_genomic_regions.txt
nuclearregionsfile=${PWD}/../scripts_and_programs/nuclear_landmarks.txt
FISH_data=${PWD}/../scripts_and_programs/FISH_in_2014_09_14_MolCell_Grob_Grossniklaus.txt
nlines=84511
pi=3.14159265359

for dir in XXXwdirXXX ;
do
    echo $dir
    cd $dir

    mkdir _tmp_IHI
    cd _tmp_IHI

    infile=${maindir}/${dir}/position_of_interesting_genomic_regions_for_${dir}.txt
    outnumber=${maindir}/${dir}/number_of_distinct_IHIs_clusters_using_R_for_${dir}.txt
    outoverlap=${maindir}/${dir}/overlapping_IHIs_using_R_for_${dir}.txt

    if [[ -s ${outnumber} ]];
    then
	echo "#Replica timestep threshold Nclusters" > ${outnumber}
    fi
    if [[ -s ${outoverlap} ]];
    then
	echo "#Replica timestep threshold Nclusters" > ${outoverlap}
    fi

    pwd

    replicas=$(awk '{print $1}' <(grep -v "#" ${infile}) | uniq | sort | uniq)
    #echo $nreplicas

    timesteps=$(awk '{print $NF}' <(grep -v "#" ${infile}) | uniq | sort | uniq)
    #echo $timesteps

    regions=$(awk '{print $2}' <(grep -v "#" ${infile}) | grep IHI | uniq | sort | uniq)
    nregions=$(echo $regions | awk '{print NF}')
    #echo $regions
    for replica in $replicas ; #$(seq 1 1 ${nreplicas});
    do
	#echo $replica
	for timestep in $timesteps;
	do
	    rtcheck=$(awk -v t=${timestep} -v r=${replica} 'BEGIN{c=0}{if($1==r && $2==t) c++}END{print c}' ${outoverlap})
	    if [[ $rtcheck -gt 0 ]];
	    then
		continue
	    fi

	    #echo $timestep
	    awk -v t=${timestep} -v r=${replica} '{if($1==r && $NF==t) print $0}' ${infile} | grep IHI > _${replica}_${timestep}_IHI
	    #cat _${replica}_${timestep}_IHI
	    
	    n=0
	    for i in $(seq 1 1 ${nregions});
	    do
		region1=$(echo $regions | awk -v i=$i '{print $i}')

		check=$(echo $region1 | grep IHI | wc -l)
		if [[ $check -eq 0 ]];
		then
		    continue
		fi
		#echo $region1
		#25 IHI7B 15.0397 70.8257 -13.6042 12.8703 19000000
		x=$( awk -v c=${region1} '{if($2==c) print $3}' _${replica}_${timestep}_IHI)
		y=$( awk -v c=${region1} '{if($2==c) print $4}' _${replica}_${timestep}_IHI)
		z=$( awk -v c=${region1} '{if($2==c) print $5}' _${replica}_${timestep}_IHI)
		Rg=$(awk -v c=${region1} '{if($2==c) print $6}' _${replica}_${timestep}_IHI)
		#echo $x $y $z $Rg

                #Replica timestep Region1 Region2 Rg1 Rg2 COMDist_12 Overlap
		for j in $(seq $i 1 ${nregions});
		do
		    region2=$(echo $regions | awk -v j=$j '{print $j}')
		    check=$(echo $region2 | grep IHI | awk -v r1=${region1} '{if($1!=r1) print $0}' | wc -l)
		    if [[ $check -eq 0 ]];
		    then
			continue
		    fi
		    n=$(($n+1))
		    #echo "Pair $region1 $region2 $n"

		    awk -v r=${replica} -v t=${timestep} -v pi=$pi -v c=$region1 -v x=$x -v y=$y -v z=$z -v Rg=$Rg '{R=Rg*sqrt(5/3); r=$6*sqrt(5/3); if(R<=r){V=4./3.*R*R*R}else{V=4./3.*r*r*r}; d=sqrt((x-$3)*(x-$3)+(y-$4)*(y-$4)+(z-$5)*(z-$5)); if(0<d && d<=(R+r)){o=(R+r-d)*(R+r-d)*(d*d + 2*d*r + 2*d*R - 3*R*R - 3*r*r + 6*r*R)/(12 * d * V)}; if(d==0){o=1.0}; if(d>(R+r)){o=0}; if(d>0.0 && o>0)print r,t,c,$2,R,r,d,o}' <( awk -v r1=${region1} -v r2=${region2} '{if($2==r2) print $0}' _${replica}_${timestep}_IHI ) >> _overlap_${replica}_${timestep}_IHI
		    #(Rg+$6-d)*(Rg+$6-d)*(d*d+2.*d*$6+2*d*Rg-3.*Rg*Rg-3.*$6*$6+6.*$6*Rg)/(12.*d*V)
		done # Close cycle over $region2
	    done # Close cycle over $region2
	    cat _overlap_${replica}_${timestep}_IHI >> ${outoverlap}

	    # Count the number of distinct IHIs as a function of the overlap threshold
            #echo $infile ${replica} ${regions} ${timesteps}

	    awk '{print $3,$4,$NF}' _overlap_${replica}_${timestep}_IHI | sed -e "s/IHI1A/0/g" -e "s/IHI1B/1/g" -e "s/IHI2A/2/g" -e "s/IHI2B/3/g" -e "s/IHI3A/4/g" -e "s/IHI3B/5/g" -e "s/IHI4A/6/g" -e "s/IHI4B/7/g" -e "s/IHI5A/8/g" -e "s/IHI5B/9/g" -e "s/IHI6A/10/g" -e "s/IHI6B/11/g" -e "s/IHI7A/12/g" -e "s/IHI7B/13/g" -e "s/IHI8A/14/g" -e "s/IHI8B/15/g" -e "s/IHI9A/16/g" -e "s/IHI9B/17/g" -e "s/IHI10A/18/g" -e "s/IHI10B/19/g" > _overlaps

	    #cat _overlaps
	    for threshold in 0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 ;
	    do
		#python ${scriptsdir}/IHIs_KEEs_mutual_distances_cCNAG.py $threshold $nregions
		echo $replica $timestep $threshold $(python ${scriptsdir}/IHIs_KEEs_mutual_distances_cCNAG.py $threshold $nregions) >> ${outnumber}
	    done
	    rm _${replica}_${timestep}_IHI _overlap_${replica}_${timestep}_IHI _overlaps

	done # Close cycle over $timestep
    done # Close cycle over $replica
    cd .. # Exit _tmp_IHI
    rm -r _tmp_IHI
    cd .. # Exit $dir    
    
done # Close cycle over directories
echo "Analysis DONE!"