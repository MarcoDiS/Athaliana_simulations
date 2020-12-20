lmpfile=1_create_mitotic_links_linear_case.lmp

for replica in $(seq 1 1 50);
do

    if [[ -d replica_${replica} ]];
    then
	continue
    fi

    mkdir -p replica_${replica}
    cd replica_${replica}
    echo replica_${replica}

    cp ../../0_generate_initial_rosettes_conformations/Initial_conformation_replica_${replica}.dat .
    cat Initial_conformation_replica_${replica}.dat | awk '{if($3=="xlo" || $3=="ylo" || $3=="zlo"){print -90.00,90.00,$3,$4; next}; print $0}' > _a ; mv _a Initial_conformation_replica_${replica}.dat
    seed=$(od -An -N3 -i /dev/urandom)

    echo "Replica ${replica} $seed" >> ../used_seed.log

    sed -e "s/XXXseedXXX/${seed}/g" -e "s/XXXreplicaXXX/${replica}/g" ../${lmpfile} > simulation_replica_${replica}.lmp

    rm -fr *XYZ
    time ( mpirun -np 1 ~/LAMMPS/lammps-31Mar17_parallel_version/src/lmp_mpi -in simulation_replica_${replica}.lmp -log none -echo screen )

    cd ..

done
