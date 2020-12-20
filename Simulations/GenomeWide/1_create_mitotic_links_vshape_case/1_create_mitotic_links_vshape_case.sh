lmpfile=1_create_mitotic_links_vshape_case.lmp

for replica in $(seq 1 1 50);
do

    if [[ -d replica_${replica} ]];
    then
	continue
    fi

    mkdir -p replica_${replica}
    cd replica_${replica}
    echo replica_${replica}

    seed=$(od -An -N3 -i /dev/urandom)

    sed -e "s/XXXseedXXX/${seed}/g" -e "s/XXXreplicaXXX/${replica}/g" ../${lmpfile} > simulation_replica_${replica}.lmp

    rm -fr *XYZ

    time ( mpirun -np 1 ~/LAMMPS/lammps-31Mar17_parallel_version/src/lmp_mpi -in simulation_replica_${replica}.lmp -log none -echo screen )

    cd ..

done
