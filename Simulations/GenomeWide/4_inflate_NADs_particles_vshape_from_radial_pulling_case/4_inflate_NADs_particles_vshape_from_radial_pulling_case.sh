lmpfile=4_inflate_NADs_particles_vshape_from_radial_pulling_case.lmp
run=500000
dump=10000
R=$(awk 'BEGIN{print 83.3333+100}')

for replica in $(seq 1 1 50);
do
    replicadir=replica_${replica}
    
    if [[ -d ${replicadir} ]];
    then
	continue
    fi
    
    if [[ ! -e ../2_radial_pulling_vshape_case/replica_${replica}/2_radial_pulling_vshape_case_part2_replica_${replica}.txt ]];
    then
	continue
    fi

    com=$(cat <(grep NOR ../0_generate_initial_rosettes_conformations/NADs_particles/NADs_particles_NORs_only.txt | awk '{print $3}' | sed "s/*/ /g") <(tail -84502 ../1_create_mitotic_links_vshape_case/replica_${replica}/chrs10_system_preparation_mitotic_links_replica_${replica}.XYZ) | awk '{if(NF==2){s1[NR]=$1;s2[NR]=$2}else{for(i in s1){if(s1[i]<=$1&&$1<=s2[i]) print $0}}}' | awk '{x+=$3;y+=$4;z+=$5;cnt++}END{print x/cnt,y/cnt,z/cnt}')
    xcom=$(echo $com | awk '{print $1}')
    ycom=$(echo $com | awk '{print $2}')
    zcom=$(echo $com | awk '{print $3}')

    echo $replicadir
    
    mkdir -p ${replicadir}
    cd ${replicadir}
    cp ../${lmpfile} simulation_replica_${replica}.lmp

    # Translate all the system by (0,0,0)    
    cat ../../2_radial_pulling_vshape_case/replica_${replica}/2_radial_pulling_vshape_case_part2_replica_${replica}.txt | awk -v x=${xcom} -v y=${ycom} -v z=${zcom} '{if(NF==9){print $1,$2,$3,$4,$5,$6,$7,$8,$9;next};{if($3=="xlo"){print x-90.00,x+90.00,$3,$4;next};if($3=="ylo"){print y-90.00,y+90.00,$3,$4;next};if($3=="zlo"){print z-90.00,z+90.00,$3,$4;next};print $0}}' > simulation_replica_${replica}    
    mv simulation_replica_${replica} 2_radial_pulling_vshape_case_part2_replica_${replica}.txt
    
    seed=$(od -An -N3 -i /dev/urandom)
    kappa=50.
    
    #32767
    #R=0.0
    sed -e "s/XXXdumpXXX/${dump}/g" -e "s,XXXrunXXX,${run},g" -e "s/XXXxNADsXXX/${xcom}/g" -e "s/XXXyNADsXXX/${ycom}/g" -e "s/XXXzNADsXXX/${zcom}/g" -e "s,XXXRXXX,${R},g" -e "s,XXXycenterXXX,${ycenter},g" -e "s,XXXkappaXXX,${kappa},g" -e "s/XXXseedXXX/${seed}/g" -e "s/XXXreplicaXXX/${replica}/g" simulation_replica_${replica}.lmp > simulation_replica_${replica}
    mv simulation_replica_${replica} simulation_replica_${replica}.lmp
    
    rm -fr *XYZ

    time ( mpirun -np 8 ~/LAMMPS/lammps-31Mar17_parallel_version/src/lmp_mpi -in simulation_replica_${replica}.lmp -log none -echo screen )
    
    cd .. # Exit ${replicadir}
    
done
