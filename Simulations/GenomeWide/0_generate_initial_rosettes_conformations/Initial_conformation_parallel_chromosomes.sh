
radius=166.6666
for replica in $(seq 1 1 50);
do

    outfile=Initial_conformation_replica_${replica}.dat
    NADsfile=NADs_particles/NADs_particles_NORs_only.txt
    inangles=NADs_particles/angles_with_NADs_particles.txt
    inbonds=NADs_particles/bonds_with_NADs_particles.txt

    nangles=$(tail -1 ${inangles}  | awk '{print $1}')
    nbonds=$(tail -1  ${inbonds} | awk '{print $1}')

    echo $seed $nangles $nbonds

    check=1
    
    while [[ $check -eq 1 ]];
    do
	seed=$(od -An -N3 -i /dev/random)

	./Initial_conformation_parallel_chromosomes -a ${radius} -b ${radius} -c ${radius} -i 15.00 -e 40.00 -n 0 -s ${seed} -m 1 > _coords 2> /dev/null #output_replica_${replica}.log

    done

    cat _coords | awk -v na=${nangles} '{if($2=="angles"){print na,$2}else{print $0}}' > _a ; mv _a _coords
    cat _coords | awk '{if($3=="ylo"){print -300,300,$3,$4}else{print $0}}' > _a ; mv _a _coords
    cat _coords | awk '{if($3=="xlo"){print -300,300,$3,$4}else{print $0}}' > _a ; mv _a _coords
    cat _coords | awk '{if($3=="zlo"){print -300,300,$3,$4}else{print $0}}' > _a ; mv _a _coords
    cat _coords ${inbonds} ${inangles} > ${outfile}
    
    rm -fr _coords

    echo "replica_${replica} ${seed}" >> used_seeds.txt

done
