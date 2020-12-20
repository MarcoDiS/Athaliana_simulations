lmpfileheader=$(ls -1 5_prodruns_linear_case_template_header.lmp | grep -v NEUT)
lmpfileinteractions=$(ls -1 5_prodruns_linear_case_template_interactions.lmp | grep -v NEUT)
lmpfiletail=$(ls -1 5_prodruns_linear_case_template_tail.lmp | grep -v NEUT)

#Parameters of the NAD interactions depend on the NADs particles radius
sigma1=0.5
sigma5=2.2 ; # 0.64 of the nuclear volume expected for Random CLose Packing (RCP)

sigma15=$(awk -v s1=${sigma1} -v s5=${sigma5} 'BEGIN{print s1+s5}')
rc15=$(awk -v s15=${sigma15} 'BEGIN{print s15*1.12246152962189}')
R015=$(awk -v s15=${sigma15} 'BEGIN{print s15*1.5}')
echo "1-5 : ${sigma15} ${rc15} ${R015}"

sigma55=$(awk -v s1=${sigma1} -v s5=${sigma5} 'BEGIN{print s5+s5}')
# The NADs are always attractive between each other                                                            
rc55=$(awk -v s55=${sigma55} 'BEGIN{print s55*2.5}')
R055=$(awk -v s55=${sigma55} 'BEGIN{print s55*1.5}')
echo "5-5 : ${sigma55} ${rc55} ${R055}"

#Parameters of the run
dump=500000
run=20000000
timestep=0.006
#######

scenarioNAD=NAD_ATTR1.0andREPU1.0
NADattraction=1.0 # They are packed together to reach the expected size of the nucleolus
sNAD="NAD_ATTR${NADattraction}andREPU${NADattraction}"

nreplicas=10

#Optimal for CH_NEUTandREPU 0.000390625
optimalCHrepulsion=0.000390625 #0.0001953125 0.000390625 0.00078125 0.0015625 0.003125 0.00625 0.0125 0.025 0.05 0.50 1.0
optimalACattraction=0.20       #0.025 0.0375 0.05 0.10 0.20 0.40
optimalFHattraction=0.50       #0.005 0.50 0.75 1.0
#Optimal for CH_ATTRandNEUT 0.30
optimalCHattraction=0.30 #0.05 0.10 0.20 0.30 0.40 0.50 1.0
#optimalACattraction=0.50 #0.025 0.0375 0.05 0.10 0.20 0.40 0.50 1.0
#optimalFHattraction=0.50       #0.005 0.50 1.0
#Optimal for CH_ATTRandREPU 0.000390625
#optimalCHrepulsion=0.00001 #0.0001953125 0.000390625 0.00078125 0.0015625 0.003125 0.00625 0.0125 0.025 0.05 0.50 1.0
#optimalCHattraction=0.20   #0.025 0.0375 0.05 0.10 0.20 0.25 0.40
#optimalACattraction=0.50       #0.005 0.50 1.0
#optimalFHattraction=0.50       #0.005 0.50 1.0


for CHattraction in ${optimalCHattraction} ; 
do
    for CHrepulsion in ${optimalCHrepulsion} ;
    do
	for ACattraction in ${optimalACattraction} ; 
	do
	    for FHattraction in 0.5 ;
	    do
		for scenarioFH in FH_ATTRandNEUT ; # FH_NEUTandNEUT FH_ATTRandNEUT
		do   	    
		    for scenarioAC in AC_ATTRandNEUT ; #AC_NEUTandREPU AC_ATTRandNEUT AC_NEUTandNEUT ;
		    do   	    
			for scenarioCH in CH_NEUTandREPU ; # CH_ATTRandNEUT CH_NEUTandNEUT CH_ATTRandREPU CH_NEUTandREPU ;
			do
		            # FH scenario
			    if [[ ${scenarioFH} == "FH_NEUTandNEUT" ]];
			    then
				sFH="FH_NEUTandNEUT"
			    fi
			    if [[ ${scenarioFH} == "FH_ATTRandNEUT" ]];
			    then
				sFH="FH_ATTR${FHattraction}andNEUT"
			    fi
   		            # AC scenario
			    if [[ ${scenarioAC} == "AC_ATTRandNEUT" ]];
			    then
				sAC="AC_ATTR${ACattraction}andNEUT"
			    fi
			    if [[ ${scenarioAC} == "AC_NEUTandNEUT" ]];
			    then
				sAC="AC_NEUTandNEUT"
			    fi
  		            # CH scenario
			    if [[ ${scenarioCH} == "CH_ATTRandREPU" ]];
			    then
				sCH="CH_ATTR${CHattraction}andREPU${CHrepulsion}"
			    fi
			    if [[ ${scenarioCH} == "CH_ATTRandNEUT" ]];
			    then
				sCH="CH_ATTR${CHattraction}andNEUT"
			    fi
			    if [[ ${scenarioCH} == "CH_NEUTandREPU" ]];
			    then
				sCH="CH_NEUTandREPU${CHrepulsion}"
			    fi
			    if [[ ${scenarioCH} == "CH_NEUTandNEUT" ]];
			    then
				sCH="CH_NEUTandNEUT"
			    fi
			    
			    wdir=system_${sNAD}_${sAC}_${sFH}_${sCH}
			    echo $wdir
			    
			    mkdir -p ${wdir}
			    cd ${wdir}
			    echo ${wdir}
			    
			    for replica in $(seq 1 1 $nreplicas);
			    do
				replicadir=replica_${replica}
				if [[ -d ${replicadir} ]];
				then
				    continue
				fi
				
				mkdir -p ${replicadir}
				cd ${replicadir}
				echo $replicadir
				
				nparticles=30975
				
				seed=$(od -An -N3 -i /dev/urandom)
				
				sed -e "s/XXXdumpXXX/${dump}/g" -e "s/XXXtimestepXXX/${timestep}/g" -e "s/XXXrunXXX/${run}/g" -e "s/XXXseedXXX/${seed}/g" -e "s/XXXsigma55XXX/${sigma55}/g" -e "s/XXXsigma15XXX/${sigma15}/g" -e "s/XXXrc55XXX/${rc55}/g" -e "s/XXXrc15XXX/${rc15}/g" -e "s/XXXR055XXX/${R055}/g" -e "s/XXXR015XXX/${R015}/g" -e "s/XXXxNADsXXX/${xcom}/g" -e "s/XXXyNADsXXX/${ycom}/g" -e "s/XXXzNADsXXX/${zcom}/g" -e "s/XXXNADattractionXXX/${NADattraction}/g" -e "s/XXXrepulsionXXX/${repulsion}/g" -e "s/XXXattractionXXX/${attraction}/g" -e "s/#${scenarioNAD}//g" -e "s/#${scenarioAC}//g" -e "s/#${scenarioFH}//g" -e "s/#${scenarioCH}//g" -e "s/XXXreplicaXXX/${replica}/g" -e "s/XXXradiusXXX/${sigma5}/g" ../../${lmpfileheader} | grep -v scenario >> _tmp.lmp

				sed -e "s/XXXdumpXXX/${dump}/g" -e "s/XXXtimestepXXX/${timestep}/g" -e "s/XXXrunXXX/${run}/g" -e "s/XXXreplicaXXX/${replica}/g" -e "s/XXXseedXXX/${seed}/g" -e "s/XXXsigma55XXX/${sigma55}/g" -e "s/XXXsigma15XXX/${sigma15}/g" -e "s/XXXrc55XXX/${rc55}/g" -e "s/XXXrc15XXX/${rc15}/g" -e "s/XXXR055XXX/${R055}/g" -e "s/XXXR015XXX/${R015}/g" -e "s/XXXxNADsXXX/${xcom}/g" -e "s/XXXyNADsXXX/${ycom}/g" -e "s/XXXzNADsXXX/${zcom}/g" -e "s/XXXNADattractionXXX/${NADattraction}/g" -e "s/XXXCHrepulsionXXX/${CHrepulsion}/g" -e "s/XXXCHattractionXXX/${CHattraction}/g" -e "s/XXXACrepulsionXXX/${ACrepulsion}/g" -e "s/XXXACattractionXXX/${ACattraction}/g" -e "s/XXXFHrepulsionXXX/${FHrepulsion}/g" -e "s/XXXFHattractionXXX/${FHattraction}/g" -e "s/#${scenarioNAD}//g" -e "s/#${scenarioAC}//g" -e "s/#${scenarioFH}//g" -e "s/#${scenarioCH}//g" -e "s/XXXradiusXXX/${sigma5}/g" ../../${lmpfileinteractions} | grep -v "#CH\|#AC\|#FH" | grep -v scenario >> _tmp.lmp

				sed -e "s/XXXdumpXXX/${dump}/g" -e "s/XXXtimestepXXX/${timestep}/g" -e "s/XXXrunXXX/${run}/g" -e "s/XXXreplicaXXX/${replica}/g" -e "s/XXXseedXXX/${seed}/g" -e "s/XXXsigma55XXX/${sigma55}/g" -e "s/XXXsigma15XXX/${sigma15}/g" -e "s/XXXrc55XXX/${rc55}/g" -e "s/XXXrc15XXX/${rc15}/g" -e "s/XXXR055XXX/${R055}/g" -e "s/XXXR015XXX/${R015}/g" -e "s/XXXxNADsXXX/${xcom}/g" -e "s/XXXyNADsXXX/${ycom}/g" -e "s/XXXzNADsXXX/${zcom}/g" -e "s/XXXNADattractionXXX/${NADattraction}/g" -e "s/XXXrepulsionXXX/${repulsion}/g" -e "s/XXXattractionXXX/${attraction}/g" -e "s/#${scenarioNAD}//g" -e "s/#${scenarioAC}//g" -e "s/#${scenarioFH}//g" -e "s/#${scenarioCH}//g" -e "s/XXXradiusXXX/${sigma5}/g" ../../${lmpfiletail} | grep -v scenario >> _tmp.lmp
			    
				rm -fr *XYZ

				time ( mpirun -np 4 ~/LAMMPS/lammps-31Mar17_parallel_version/src/lmp_mpi -in _tmp.lmp -log none -echo screen )
				
				cd .. # Exit ${replicadir}
				
			    done
			    cd .. ${wdir} # Exit ${wdir}
			done # Close cycle over ${CHscenario}
		    done # Close cycle over ${ACscenario}
		done # Close cycle over ${FHscenario}
	    done # Close cycle over ${FHattraction}
	done # Close cycle over ${ACattraction}
    done # Close cycle over ${CHrepulsion}
done # Close cycle over ${CHattraction}