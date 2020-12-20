ndraws=50
rm -fr epigenomic_state_of_IHIs.txt _tmp epigenomic_state_of_randomised_IHIs.txt
for IHI in $(awk '{print $1}' IHIs_positions_with_chromosome.txt);
do
    chr=$(awk -v IHI=${IHI} '{if($1==IHI) print $4}' IHIs_positions_with_chromosome.txt)
    s1=$(awk  -v IHI=${IHI} '{if($1==IHI) print $2}' IHIs_positions_with_chromosome.txt)
    s2=$(awk  -v IHI=${IHI} '{if($1==IHI) print $3}' IHIs_positions_with_chromosome.txt)
    l=$(($s2-$s1))
    chrLength=$(awk -v c=$chr -v l=$l '{if($1==c) print $NF-l}' INPUT_PARAMETERS_COMPLETE.txt)
    echo $IHI $chr $chrLength $s1 $s2 $l

    # IHI
    awk -v s1=$s1 -v s2=$s2 '{if(s1<=$1 && $1<=s2) print $2}' chromatin_states_genome_wide.txt >> epigenomic_state_of_IHIs.txt 
    nl=$(awk -v s1=$s1 -v s2=$s2 '{if(s1<=$1 && $1<=s2) print $2}' chromatin_states_genome_wide.txt | wc -l | awk '{print $1}')

    # Random regions
    offset=$(awk 'BEGIN{s=0}{print $1,s; s+=$NF}' INPUT_PARAMETERS_COMPLETE.txt | awk -v c=$chr '{if($1==c) print $2}')
    for i in $(seq 1 1 $ndraws);
    do
	s1=$(( $offset + ( RANDOM % $chrLength )  + 1 ))
	s2=$(($l+$s1))
	echo $IHI $chr $offset $chrLength $s1 $s2 $l 
	awk -v s1=$s1 -v s2=$s2 -v i=$i '{if(s1<=$1 && $1<=s2) print i,$2}' chromatin_states_genome_wide.txt > _tmp 
	nlines=$(wc -l _tmp | awk '{print $1}')
	if [[ $nlines -ne $nl ]];
	then
	    cat _tmp
	    echo $nlines $nl
	    exit
	fi
	cat _tmp >> epigenomic_state_of_randomised_IHIs.txt 
    done # Close cycle over $ndraws
    wc -l epigenomic_state_of_randomised_IHIs.txt
done # Close cycle over $IHI

rm -fr _tmp
awk '{h[$1]++}END{for(i in h) print i,h[i]}' epigenomic_state_of_IHIs.txt > histrogram_of_epigenetic_states_of_IHIs.txt
for i in $(seq 1 1 $ndraws);
do
    awk -v i=${i} '{if($1==i)h[$2]++}END{for(i in h) print i,h[i]}' epigenomic_state_of_randomised_IHIs.txt >> _tmp
done # Close cycle over $ndraws
awk -v nd=$ndraws '{h[$1]+=$2; h2[$1]+=$2*$2; cnt[$1]++}END{for(i in h){avg=h[i]/nd;avg2=h2[i]/nd;stddev=sqrt(avg2-avg*avg); print i,avg,stddev,nd}}' _tmp > histrogram_of_epigenetic_states_of_randomised_IHIs.txt