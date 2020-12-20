#awk -v nc=3 'BEGIN{print (10143*3000/(30*30*30)/0.00441745*nc)^(1/3); print nc*10143}'

from pytadbit.modelling.lammps_modelling import generate_chromosome_rosettes_conformation_with_pbc

# From 52795 to 58989 both included
nparticles = (7529-1334)
nchrs = 5
chromosome_particle_numbers=[nparticles]*nchrs

side = (nparticles*3000*nchrs/(30*30*30)/0.00441745)**(1./3.)

print "%d Chromosomes of length %d in a cube of side %s" % (nchrs, nparticles, side)

r = XXXreplicaXXX
seed = XXXseedXXX

for replica in xrange(r, r+1,1):
    initial_conformation = "Initial_rosette_conformation_with_pbc_replica_%s.dat" % (replica)
    generate_chromosome_rosettes_conformation_with_pbc ( chromosome_particle_numbers, fractional_radial_positions=None,
                                                         confining_environment=['cube',side], rosette_radius=12.0 , particle_radius=0.5 ,
                                                         seed_of_the_random_number_generator=seed ,
                                                         number_of_conformations=1,
                                                         outfile = initial_conformation)
