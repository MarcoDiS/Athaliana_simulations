IFS=$'\n'  

#INPUT
#chr1A 30427671 10143
#chr1B 30427671 10143

#OUTPUT
#CHR   PART  REMAP
#chr1A 1     1
#
#chr5B 79432 39716

awk 'BEGIN{offset=0}{for(i=1;i<=$3;i++){print $1"A",i+offset,i+offset/2;print $1"B",i+offset+$3,i+offset/2;};offset+=$3*2}' <( grep -v "#" INPUT_PARAMETERS_COMPLETE.txt | grep A | sed "s/A//g") | sort -k 2n
