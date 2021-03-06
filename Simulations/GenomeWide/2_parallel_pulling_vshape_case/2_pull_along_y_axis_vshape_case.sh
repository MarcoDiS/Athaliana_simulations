lmpfile=2_pull_along_y_axis_vshape_case.lmp

for replica in $(seq 1 1 50);
do
    
    if [[ -d replica_${replica} ]];
    then
	continue
    fi
    
    if [[ ! -e ../1_create_mitotic_links_vshape_case/replica_${replica}/1_create_mitotic_links_replica_${replica}.txt ]];
    then
	continue
    fi
    replicadir=replica_${replica}
    echo $replicadir
    
    mkdir -p ${replicadir}
    cd ${replicadir}
    cp ../${lmpfile} simulation_replica_${replica}.lmp
    for c in $(cat ../centromeres.txt | awk '{print $(NF-1)"_"$NF"_"$2}');
    do
	s1=$(echo  $c | sed -e "s/_/ /g" | awk '{print $1}')
	s2=$(echo  $c | sed -e "s/_/ /g" | awk '{print $2}')
	cen=$(echo $c | sed -e "s/_/ /g" | awk '{print $3}')
	echo $c $s1 $s2 $cen
	x=$(awk -v s1=$s1 -v s2=$s2 '{if(NF==9 && s1<=$1 && $1<=s2){x+=$4;c++}}END{print x/c}' ../../1_create_mitotic_links_vshape_case/replica_${replica}/1_create_mitotic_links_replica_${replica}.txt)
	y=$(awk -v s1=$s1 -v s2=$s2 '{if(NF==9 && s1<=$1 && $1<=s2){x+=$5;c++}}END{print x/c}' ../../1_create_mitotic_links_vshape_case/replica_${replica}/1_create_mitotic_links_replica_${replica}.txt)
	z=$(awk -v s1=$s1 -v s2=$s2 '{if(NF==9 && s1<=$1 && $1<=s2){x+=$6;c++}}END{print x/c}' ../../1_create_mitotic_links_vshape_case/replica_${replica}/1_create_mitotic_links_replica_${replica}.txt)
	echo $x $y $z

	sed -e "s/XXXx${cen}XXX/${x}/g" -e "s/XXXy${cen}XXX/${y}/g" -e "s/XXXz${cen}XXX/${z}/g" simulation_replica_${replica}.lmp > simulation_replica_${replica}
	mv simulation_replica_${replica} simulation_replica_${replica}.lmp
    done

    # Translate all the system by (0,+150,0) in the y direction
    cat ../../1_create_mitotic_links_vshape_case/replica_${replica}/1_create_mitotic_links_replica_${replica}.txt | awk '{if(NF==9){print $1,$2,$3,$4,$5+150,$6,$7,$8,$9;next};{if($3=="xlo" || $3=="ylo" || $3=="zlo"){print -300.00,300.00,$3,$4;next};print $0}}' > simulation_replica_${replica}    
    mv simulation_replica_${replica} 1_create_mitotic_links_replica_${replica}.txt
    
    seed=$(od -An -N3 -i /dev/urandom)
    kappa=50.
    
    R=0.0
    sed -e "s,XXXRXXX,${R},g" -e "s,XXXycenterXXX,${ycenter},g" -e "s,XXXkappaXXX,${kappa},g" -e "s/XXXseedXXX/${seed}/g" -e "s/XXXreplicaXXX/${replica}/g" simulation_replica_${replica}.lmp > simulation_replica_${replica}
    mv simulation_replica_${replica} simulation_replica_${replica}.lmp
    
    rm -fr *XYZ

    time ( mpirun -np 1 ~/LAMMPS/lammps-31Mar17_parallel_version/src/lmp_mpi -in simulation_replica_${replica}.lmp -log none -echo screen )
    
    cd .. # Exit ${replicadir}
    
done
