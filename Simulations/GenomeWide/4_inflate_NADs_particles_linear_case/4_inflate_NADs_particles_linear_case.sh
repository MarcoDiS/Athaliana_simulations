r=$1
radius=$2

lmpfile=4_inflate_NADs_particles_linear_case.lmp
echo $lmpfile

for replica in $r ;
do

    replicadir=replica_${replica}_radius_${radius}
    
    mkdir -p ${replicadir}
    cd ${replicadir}
    echo $replicadir
    
    rm -fr *XYZ
    sigma1=0.5
    
    seed=$(od -An -N3 -i /dev/random)
    echo $seed
	
    sigma4=0.5
    
    sigma14=$(awk -v s1=${sigma1} -v s4=${sigma4} 'BEGIN{print s1+s4}')
    rc14=$(awk -v s14=${sigma14} 'BEGIN{print s14*1.12246152962189}')
    R014=$(awk -v s14=${sigma14} 'BEGIN{print s14*1.5}')
    echo "1-4 : ${sigma14} ${rc14} ${R014}"
    
    sigma44=$(awk -v s1=${sigma1} -v s4=${sigma4} 'BEGIN{print s4+s4}')
    rc44=$(awk -v s44=${sigma44} 'BEGIN{print s44*2.5}')
    R044=$(awk -v s44=${sigma44} 'BEGIN{print s44*1.5}')
    echo "4-4 : ${sigma44} ${rc44} ${R044}"
    
    run=1000000
    dump=10000
    timestep=0.006
    sed -e "s/XXXdumpXXX/${dump}/g" -e "s/XXXtimestepXXX/${timestep}/g" -e "s/XXXrunXXX/${run}/g" -e "s/XXXreplicaXXX/${replica}/g" -e "s/XXXseedXXX/${seed}/g" -e "s/XXXsigma44XXX/${sigma44}/g" -e "s/XXXsigma14XXX/${sigma14}/g" -e "s/XXXrc44XXX/${rc44}/g" -e "s/XXXrc14XXX/${rc14}/g" -e "s/XXXR044XXX/${R044}/g" -e "s/XXXR014XXX/${R014}/g" -e "s/XXXxNADsXXX/${xcom}/g" -e "s/XXXyNADsXXX/${ycom}/g" -e "s/XXXzNADsXXX/${zcom}/g" ../$lmpfile > simulation_replica_${replica}.lmp
    
    rm -fr output_lammps_replica_1.log
    mpirun -np 8 ~/LAMMPS/lammps-31Mar17_parallel_version/src/lmp_mpi -in simulation_replica_${replica}.lmp -log none -echo screen > output_lammps_replica_1.log
    
    for sigma4 in $(seq 0.5 0.1 ${radius});
    do

	echo $replica $sigma4 "Rg: "$(bash ~/compute_Rg.sh)
	
	sigma14=$(awk -v s1=${sigma1} -v s4=${sigma4} 'BEGIN{print s1+s4}')
	rc14=$(awk -v s14=${sigma14} 'BEGIN{print s14*1.12246152962189}')
	R014=$(awk -v s14=${sigma14} 'BEGIN{print s14*1.5}')
	echo "1-4 : ${sigma14} ${rc14} ${R014}"
	
	sigma44=$(awk -v s1=${sigma1} -v s4=${sigma4} 'BEGIN{print s4+s4}')
	#rc44=$(awk -v s44=${sigma44} 'BEGIN{print s44*1.12246152962189}')
	rc44=$(awk -v s44=${sigma44} 'BEGIN{print s44*2.5}')
	R044=$(awk -v s44=${sigma44} 'BEGIN{print s44*1.5}')
	echo "4-4 : ${sigma44} ${rc44} ${R044}"
	
	run=1000
	dump=1000
	timestep=0.0001
	sed -e "s/XXXdumpXXX/${dump}/g" -e "s/XXXtimestepXXX/${timestep}/g" -e "s/XXXrunXXX/${run}/g" -e "s/XXXreplicaXXX/${replica}/g" -e "s/#read_restart/read_restart/g" -e "s/read_data/#read_data/g" -e "s/XXXseedXXX/${seed}/g" -e "s/XXXsigma44XXX/${sigma44}/g" -e "s/XXXsigma14XXX/${sigma14}/g" -e "s/XXXrc44XXX/${rc44}/g" -e "s/XXXrc14XXX/${rc14}/g" -e "s/XXXR044XXX/${R044}/g" -e "s/XXXR014XXX/${R014}/g" -e "s/reset_timestep/#reset_timestep/g" -e "s/XXXxNADsXXX/${xcom}/g" -e "s/XXXyNADsXXX/${ycom}/g" -e "s/XXXzNADsXXX/${zcom}/g" ../${lmpfile} > simulation_replica_${replica}.lmp
	
	mpirun -np 8 ~/LAMMPS/lammps-31Mar17_parallel_version/src/lmp_mpi -in simulation_replica_${replica}.lmp -log none -echo screen >> output_lammps_replica_1.log
	
    done

    run=1000000
    dump=10000
    timestep=0.003
    sed -e "s/XXXdumpXXX/${dump}/g" -e "s/XXXtimestepXXX/${timestep}/g" -e "s/XXXrunXXX/${run}/g" -e "s/XXXreplicaXXX/${replica}/g" -e "s/#read_restart/read_restart/g" -e "s/read_data/#read_data/g" -e "s/XXXseedXXX/${seed}/g" -e "s/XXXsigma44XXX/${sigma44}/g" -e "s/XXXsigma14XXX/${sigma14}/g" -e "s/XXXrc44XXX/${rc44}/g" -e "s/XXXrc14XXX/${rc14}/g" -e "s/XXXR044XXX/${R044}/g" -e "s/XXXR014XXX/${R014}/g" -e "s/reset_timestep/#reset_timestep/g" -e "s/XXXxNADsXXX/${xcom}/g" -e "s/XXXyNADsXXX/${ycom}/g" -e "s/XXXzNADsXXX/${zcom}/g" ../${lmpfile} > simulation_replica_${replica}.lmp
	
    mpirun -np 8 ~/LAMMPS/lammps-31Mar17_parallel_version/src/lmp_mpi -in simulation_replica_${replica}.lmp -log none -echo screen >> output_lammps_replica_1.log

    cd .. # Exit ${replicadir}
done
