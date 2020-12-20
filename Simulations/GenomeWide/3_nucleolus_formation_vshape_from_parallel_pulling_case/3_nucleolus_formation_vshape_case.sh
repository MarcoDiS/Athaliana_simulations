lmpfile=3_nucleolus_formation_vshape_case.lmp
echo $lmpfile

mkdir -p ${scenario}
echo $scenario
cd ${scenario}

for replica in $(seq 1 1 50);
do
    if [[ -d replica_$replica ]];
    then
	continue
    fi
    if [[ ! -e ../../2_yaxis_pulling_vshape_case/replica_${replica}/chrs10_system_preparation_centromere_pulling_replica_${replica}.XYZ ]];
    then
	ls -lrth ../../2_yaxis_pulling_vshape_case/replica_${replica}/
	echo "I'm here"
	continue
    fi
    replicadir=replica_${replica}
    echo $replicadir
    
    com=$(cat <(grep NOR ../../0_generate_initial_rosettes_conformations/NADs_particles/NADs_particles_NORs_only.txt | awk '{print $3}' | sed "s/*/ /g") <(tail -84502 ../../2_yaxis_pulling_vshape_case/replica_${replica}/chrs10_system_preparation_centromere_pulling_replica_${replica}.XYZ) | awk '{if(NF==2){s1[NR]=$1;s2[NR]=$2}else{for(i in s1){if(s1[i]<=$1&&$1<=s2[i]) print $0}}}' | awk '{x+=$3;y+=$4;z+=$5;cnt++}END{print x/cnt,y/cnt,z/cnt}')
    xcom=$(echo $com | awk '{print $1}')
    ycom=$(echo $com | awk '{print $2}')
    zcom=$(echo $com | awk '{print $3}')
    
    mkdir -p ${replicadir}
    cd ${replicadir}
    
    seed=$(od -An -N3 -i /dev/urandom)
    sed -e "s/XXXxNADsXXX/${xcom}/g" -e "s/XXXyNADsXXX/${ycom}/g" -e "s/XXXzNADsXXX/${zcom}/g" -e "s/XXXseedXXX/${seed}/g" -e "s/XXXreplicaXXX/${replica}/g" ../../${lmpfile} > _tmp.lmp
    
    rm -fr *XYZ
    
    time ( mpirun -np 1 ~/LAMMPS/lammps-31Mar17_parallel_version/src/lmp_mpi -in _tmp.lmp -log none -echo screen )
    
    cd .. # Exit ${replicadir}
    
done
