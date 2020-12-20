# Generate the initial rod-like conformation using TADdyn
for replica in $(seq 1 1 10);
do

    seed=$(od -An -N3 -i /dev/random | awk '{print $1}')

    echo $replica $seed >> used_seeds.txt
    
    sed -e "s/XXXseedXXX/${seed}/g" -e "s/XXXreplicaXXX/${replica}/g" Initial_conformation_parallel_chromosomes.py | python
    
done
