#!/bin/bash

maindir=${PWD}
scriptsdir=${PWD}/../Analysis/
genomicregionsfile=${scriptsdir}/Utilities/interesting_genomic_regions.txt
nlines=84511

for dir in XXXwdirXXX ;
do
    echo $dir
    cd $dir
    pwd

    outfile=${maindir}/${dir}/position_of_interesting_genomic_regions_for_${dir}.txt
    touch ${outfile}    
    if [[ ! -s ${outfile} ]];
    then
	echo "#Replica GenRegion xGR yGR zGR Rg Time" > ${outfile}
    fi

    nl=$(wc -l ${outfile} | awk '{print $1}')    
    if [[ ${nl} -ne 148201 ]];
    then
	for replica in $(seq 1 1 50);
	do
	    check=$(awk -v r=$replica '{if($1==r) h[$1]++}END{for(i in h) print h[i]}' ${outfile})
	    if [[ $check -eq 2964 ]];
	    then
		echo "Replica $replica DONE!"
		continue
	    fi
	    echo replica_${replica}
	    cd contact_map/input_files/ 
	    rm -fr ~/.*~ *~ ../*~ ./*/*~ ../*/*~ ./*/*/*~ ./*/*/*/*~ ./*/*/*/*/*~ ~/*~ ~/*/*~
	    
	    for genomic_region in $(grep -v "#" ${genomicregionsfile} | awk '{print $1}' | uniq) ;
	    do
		echo ${genomic_region}
		
		s1=$(awk -v g=${genomic_region} '{if($1==g){print $2}}' ${genomicregionsfile})
		s2=$(awk -v g=${genomic_region} '{if($1==g){print $3}}' ${genomicregionsfile})
		awk -v g=${genomic_region} '{if($1==g){print $0}}' ${genomicregionsfile}
		echo $s1 $s2
		
		for file in $(ls -1 | grep XYZ | grep replica_${replica}_);
		do
		    
		    echo $file
		    
		    # Data file
		    infile=data_file.txt
    		    #echo $R $s1 $s2 $nlines $x $y $z
		    #awk -v r=${R} -v s1=${s1} -v s2=${s2} -v nl=${nlines} -v x=${x} -v y=${y} -v z=${z} '{if(NR%nl==2 && NR==2){t=$1;xc=0;yc=0;zc=0;next}; if(NR%nl==2 && NR!=2){xc=xc/cnt;yc=yc/cnt;zc=zc/cnt;dx=x-xc;dy=y-yc;dz=z-zc;d=sqrt(dx*dx+dy*dy+dz*dz); print t,xc,yc,zc,d,d/r;t=$1;xc=0;yc=0;zc=0;cnt=0; next}; if($1>=s1 && $1<=s2){xc+=$3;yc+=$4;zc+=$5;cnt++}}END{xc=xc/cnt;yc=yc/cnt;zc=zc/cnt;dx=x-xc;dy=y-yc;dz=z-zc;d=sqrt(dx*dx+dy*dy+dz*dz); print t,xc,yc,zc,d,d/r}' ${file} > ${infile}
		    awk -v s1=${s1} -v s2=${s2} -v nl=${nlines} '{if(NR%nl==2 && NR==2){t=$1;xc=0;yc=0;zc=0;next}; if(NR%nl==2 && NR!=2){xc=xc/cnt;yc=yc/cnt;zc=zc/cnt; print t,xc,yc,zc;t=$1;xc=0;yc=0;zc=0;cnt=0; next}; if($1>=s1 && $1<=s2){xc+=$3;yc+=$4;zc+=$5;cnt++}}END{xc=xc/cnt;yc=yc/cnt;zc=zc/cnt; print t,xc,yc,zc}' ${file} > ${infile}
		    
		    
		    #cat ${infile}
		    for t in $(cat ${infile} | awk '{print $1}');
		    do		    
			check=$(awk -v r=${replica} -v g=${genomic_region} -v t=${t} 'BEGIN{c=0}{if($1==r && $2==g && t==$7) c=1}END{print c}' ${outfile})
			if [[ $check -eq 1 ]];
			then
		    	    #cat ${outfile}
			    continue
			fi
			pwd
			
		        # Rg
			xc=$(awk -v t=${t} '{if($1==t) print $2}' ${infile} | uniq)
			yc=$(awk -v t=${t} '{if($1==t) print $3}' ${infile} | uniq)
			zc=$(awk -v t=${t} '{if($1==t) print $4}' ${infile} | uniq)
		        #echo ${t} ${xc} ${yc} ${zc}
			
			rg=$(awk -v t=${t} -v s1=${s1} -v s2=${s2} -v nl=${nlines} -v xc=${xc} -v yc=${yc} -v zc=${zc} 'BEGIN{f=0}{if(NR%nl==2 && t==$1){cnt=0;rg=0;rg2=0;f=1;next}; if(NR%nl==2 && t!=$1){f=0}; if($1>=s1 && $1<=s2){if(f==1){dx=xc-$3;dy=yc-$4;dz=zc-$5;r=dx*dx+dy*dy+dz*dz;rg+=r;rg2+=r*r;cnt++}}}END{avg=rg/cnt;avg2=rg2/cnt;stddev=sqrt(avg2-avg*avg);print sqrt(avg)}' ${file})
			
		        #grep -w ${t} ${infile} #| grep -w ${t}
			
			echo $replica $genomic_region $xc $yc $zc $rg $t >> ${outfile}
		        #cat ${outfile}
		    done # Close cycle over time
		    #cat ${outfile}		
		    rm -fr data_file.txt
		done # Close cycle over XYZ file
	    done # Close cycle over genomic_region
	done  # Close cycle over replicas
	cd ../../ # Exit contact_map/input_files/
	pwd
    fi
        
done # Close cycle over directories
echo "Analysis DONE!"
