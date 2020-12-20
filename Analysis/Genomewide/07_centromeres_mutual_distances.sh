#!/bin/bash

maindir=${PWD}
scriptsdir=${PWD}/../

pi=3.14159265359

for dir in XXXwdirXXX ;
do
    echo $dir
    cd $dir

    mkdir _tmp_centromere
    cd _tmp_centromere

    infile=${maindir}/${dir}/position_of_interesting_genomic_regions_for_${dir}.txt
    outnumber=${maindir}/${dir}/number_of_distinct_centromeres_clusters_using_Rg_for_${dir}.txt
    outoverlap=${maindir}/${dir}/overlapping_centromeres_using_Rg_for_${dir}.txt

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

    regions=$(awk '{print $2}' <(grep -v "#" ${infile}) | grep CEN | uniq | sort | uniq)
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
	    awk -v t=${timestep} -v r=${replica} '{if($1==r && $NF==t) print $0}' ${infile} | grep CEN > _${replica}_${timestep}
	    #cat _${replica}_${timestep}
	    
	    n=0
	    for i in $(seq 1 1 ${nregions});
	    do
		region1=$(echo $regions | awk -v i=$i '{print $i}')

		check=$(echo $region1 | grep CEN | wc -l)
		if [[ $check -eq 0 ]];
		then
		    continue
		fi
		#echo $region1
		#25 IHI7B 15.0397 70.8257 -13.6042 12.8703 19000000
		x=$( awk -v c=${region1} '{if($2==c) print $3}' _${replica}_${timestep})
		y=$( awk -v c=${region1} '{if($2==c) print $4}' _${replica}_${timestep})
		z=$( awk -v c=${region1} '{if($2==c) print $5}' _${replica}_${timestep})
		Rg=$(awk -v c=${region1} '{if($2==c) print $6}' _${replica}_${timestep})
		#echo $x $y $z $Rg

                #Replica timestep Region1 Region2 Rg1 Rg2 COMDist_12 Overlap
		for j in $(seq $i 1 ${nregions});
		do
		    region2=$(echo $regions | awk -v j=$j '{print $j}')
		    check=$(echo $region2 | grep CEN | awk -v r1=${region1} '{if($1!=r1) print $0}' | wc -l)
		    if [[ $check -eq 0 ]];
		    then
			continue
		    fi
		    n=$(($n+1))
		    #echo "Pair $region1 $region2 $n"

		    awk -v rep=${replica} -v t=${timestep} -v pi=$pi -v c=$region1 -v x=$x -v y=$y -v z=$z -v Rg=$Rg '{R=Rg; r=$6; if(R<=r){V=4./3.*R*R*R}else{V=4./3.*r*r*r}; d=sqrt((x-$3)*(x-$3)+(y-$4)*(y-$4)+(z-$5)*(z-$5)); if(0<d && d<=(R+r)){o=(R+r-d)*(R+r-d)*(d*d + 2*d*r + 2*d*R - 3*R*R - 3*r*r + 6*r*R)/(12 * d * V)}; if(d==0){o=1.0}; if(d>(R+r)){o=0}; if(d>0.0 && o>0)print rep,t,c,$2,R,r,d,o}' <( awk -v r1=${region1} -v r2=${region2} '{if($2==r2) print $0}' _${replica}_${timestep} ) >> _overlap_${replica}_${timestep}
		done # Close cycle over $region2
	    done # Close cycle over $region2
	    cat _overlap_${replica}_${timestep} >> ${outoverlap}

	    # Count the number of distinct centromeres as a function of the overlap threshold
	    awk '{print $3,$4,$NF}' _overlap_${replica}_${timestep} | sed -e "s/CEN1A/0/g" -e "s/CEN1B/1/g" -e "s/CEN2A/2/g" -e "s/CEN2B/3/g" -e "s/CEN3A/4/g" -e "s/CEN3B/5/g" -e "s/CEN4A/6/g" -e "s/CEN4B/7/g" -e "s/CEN5A/8/g" -e "s/CEN5B/9/g" > _overlaps
	    #cat _overlaps
	    for threshold in 0.00 0.10 0.20 0.30 0.31 0.32 0.33 0.34 0.35 0.36 0.37 0.38 0.39 0.40 0.50 0.60 0.70 0.80 0.90 ;
	    do
		echo $replica $timestep $threshold $(python ${scriptsdir}/Utilities/IHIs_KEEs_mutual_distances.py $threshold $nregions) >> ${outnumber}
	    done
	    rm _${replica}_${timestep} _overlap_${replica}_${timestep} _overlaps

	done # Close cycle over $timestep
    done # Close cycle over $replica
    cd .. # Exit _tmp_centromere
    rm _tmp_centromere
    cd .. # Exit $dir    
    
done # Close cycle over directories
echo "Analysis DONE!"
