OutfilePerBin=${PWD}/pvalues_occupancy_Wilcoxon_per_bin.txt
refnreplicas=50
nreplicas=50
pi=3.14159265359

dirs=$(ls -1 | grep system_N | grep -v log | grep -v cmd | grep -v IHI)

for refdir in system_NAD_ATTR1.0andREPU1.0_AC_ATTR0.20andNEUT_FH_ATTR0.50andNEUT_CH_NEUTandREPU0.000390625 system_NAD_ATTR1.0andREPU1.0_AC_NEUTandNEUT_FH_NEUTandNEUT_CH_NEUTandNEUT ;
do
    echo $refdir
    for dir in ${dirs} ;
    do
	tdir=$(echo $dir | sed "s/system_//g")
	trefdir=$(echo $refdir | sed "s/system_//g")
	redrefdir=$(echo $refdir | sed -e "s/system_//g" -e "s/NEUT/N/g" -e "s/ATTR/A/g" -e "s/REPU/R/g")
	reddir=$(echo $dir | sed -e "s/system_//g" -e "s/NEUT/N/g" -e "s/ATTR/A/g" -e "s/REPU/R/g")

	echo $dir

	if [[ ! -d ${dir}/contact_map/input_files ]];
	then
	    continue
	fi

	# Histrogram with the shells of the same volume
	for nbins in 3 5 ;
	do
	    continue
	    for state in NORs Heterochromatin Polycomb_like Active Undetermined Telomere ;
	    do
		check=$(awk -v nb=${nbins} -v state=$state -v td=${tdir} -v tr=${trefdir} '{if($1==td && $2==nb && $3==state && $NF==tr) print $0}' ${OutfilePerBin} | wc -l)
		echo $check
		if [[ $check -eq ${nbins} ]];
		then
		    echo "${dir} vs. ${refdir} ${state} DONE!"
		    continue
		fi
		
		echo $state ${nbins}

  	        # Wilcoxon test to evaluate the different mean per bin
	        # Signal data
		awk -v nr=${nreplicas} -v nbins=${nbins} -v R=2500 'BEGIN{for(i=0;i<=(nbins-1);i++){s1[i]=R*((i/nbins)^(1./3.)); s2[i]=R*(((i+1)/nbins)^(1./3.))}}{r=$1*30; for(i=0;i<=(nbins-1);i++){if(s1[i]<=r && r<=s2[i]){bin=i}}; h[$(NF-1),bin]++; cnt[$(NF-1)]++;}END{for(i=1;i<=nr;i++){for(j=0;j<=(nbins-1);j++){if(cnt[i]==0){continue}else{print s1[j],h[i,j]/cnt[i],h[i,j],cnt[i],i}}}}' ${dir}/contact_map/input_files/radial_position_${state}.txt | sort -k 1n > _tmp_signal_all_${refdir}_${dir}

	        # Control data
		awk -v nr=${refnreplicas} -v nbins=${nbins} -v R=2500 'BEGIN{for(i=0;i<=(nbins-1);i++){s1[i]=R*((i/nbins)^(1./3.)); s2[i]=R*(((i+1)/nbins)^(1./3.))}}{r=$1*30; for(i=0;i<=(nbins-1);i++){if(s1[i]<=r && r<=s2[i]){bin=i}}; h[$(NF-1),bin]++; cnt[$(NF-1)]++;}END{for(i=1;i<=nr;i++){for(j=0;j<=(nbins-1);j++){if(cnt[i]==0){continue}else{print s1[j],h[i,j]/cnt[i],h[i,j],cnt[i],i}}}}' ${refdir}/contact_map/input_files/radial_position_${state}.txt | sort -k 1n > _tmp_control_all_${refdir}_${dir}

		for bin in $(awk '{print $1}' _tmp_signal_all_${refdir}_${dir} | sort | uniq) ;
		do 
		    
		    # Signal data
		    awk -v b=${bin} '{if($1==b) print "signal", $2}'  _tmp_signal_all_${refdir}_${dir} >  _tmp_${refdir}_${dir}

		    # Control data
		    awk -v b=${bin} '{if($1==b) print "control",$2}' _tmp_control_all_${refdir}_${dir} >> _tmp_${refdir}_${dir}

		    # Entire distribution
		    sed -e "s/XXXrefdirXXX/${refdir}/g" -e "s/XXXdirXXX/${dir}/g" ${scriptsdir}/Utilities/Wilcoxon_per_bin_statistical_analysis.R > _tmp_${refdir}_${dir}.R
		    stats=$(Rscript _tmp_${refdir}_${dir}.R | grep "\[1\]" | sed "s/\"//g" | awk '{print $11}')
		    echo $nbins $bin $stats

		    mean=$(awk -v b=$bin '{if($1==b) print $2,$3}' ${dir}/contact_map/input_files/${state}_histogram_${nbins}_equal_volume_bins_radial_position.txt)
		    echo $mean
		    echo ${dir} ${nbins} ${state} ${bin} ${mean} ${stats} ${refdir} | sed "s/system_//g" >> ${OutfilePerBin}

		    #rm -fr _tmp_${refdir}_${dir}
		done # Close cycle over $bin
	        rm -fr _tmp_signal_all_${refdir}_${dir} _tmp_control_all_${refdir}_${dir} _tmp_all_${refdir}_${dir}

	    done # Close cycle over $state
	done # Close cycle over $nbins

	# Histrogram with the shells of the same thickness
	for binwidth in 1500 ; #125 500 250 ;
	do

	    for state in NORs Polycomb_like Telomere Undetermined Active Heterochromatin 
	    do
		check=$(awk -v bw=${binwidth} -v state=$state -v td=${tdir} -v tr=${trefdir} '{if($1==td && $2==bw && $3==state && $NF==tr) print $0}' ${OutfilePerBin} | wc -l)
		ntot=6
		echo $check
		if [[ $check -eq ${ntot} ]];
		then
		    echo "${dir} vs. ${refdir} ${state} DONE!"
		    continue
		fi

		echo $state ${binwidth}
		
  	        # Wilcoxon test to evaluate the different mean per bin
	        # Signal data
		awk -v pi=$pi -v nr=${nreplicas}    -v bw=${binwidth} '{bin=int(($1*30)/bw); h[$(NF-1),bin]++; cnt[$(NF-1)]++;}END{for(i=1;i<=nr;i++){for(j=0;j<=((2500-bw)/bw);j++){R=(j+1)*bw; r=(j+0)*bw; V=4./3.*pi*(R*R*R-r*r*r); if(cnt[i]==0){continue}else{print j*bw,h[i,j]/(cnt[i]),i}}}}' ${dir}/contact_map/input_files/radial_position_${state}.txt    > _tmp_signal_all

	        # Control data
		awk -v pi=$pi -v nr=${refnreplicas} -v bw=${binwidth} '{if(NF==4){bin=int(($1*30)/bw); h[$(NF-1),bin]++; cnt[$(NF-1)]++;}else{bin=int(($1*30)/bw); h[$NF,bin]++; cnt[$NF]++;}}END{for(i=1;i<=nr;i++){for(j=0;j<=((2500-bw)/bw);j++){R=(j+1)*bw; r=(j+0)*bw; V=4./3.*pi*(R*R*R-r*r*r); if(cnt[i]==0){continue}else{print j*bw,h[i,j]/(cnt[i]),i}}}}' ${refdir}/contact_map/input_files/radial_position_${state}.txt > _tmp_control_all  

		for bin in $(seq 0 ${binwidth} 2499) ;
		do 
		    
		    # Signal data
		    awk -v b=${bin} '{if($1==b) print "signal", $2}'  _tmp_signal_all_${refdir}_${dir} >  _tmp_${refdir}_${dir}

		    # Control data
		    awk -v b=${bin} '{if($1==b) print "control",$2}' _tmp_control_all_${refdir}_${dir} >> _tmp_${refdir}_${dir}
		    sed -e "s/XXXrefdirXXX/${refdir}/g" -e "s/XXXdirXXX/${dir}/g" ../scripts_and_programs/Utilities/Wilcoxon_per_bin_statistical_analysis.R > _tmp_${refdir}_${dir}.R		    
		    stats=$(Rscript _tmp_${refdir}_${dir}.R | grep "\[1\]" | sed "s/\"//g" | awk '{print $11}')
		    echo $binwidth $bin $stats

		    mean=$(awk -v b=$bin '{if($1==b) print $2,$3}' ${dir}/contact_map/input_files/${state}_histogram_at_${binwidth}nm_radial_position.txt)
		
		    echo ${dir} ${binwidth} ${state} ${bin} ${mean} ${stats} ${refdir} | sed "s/system_//g" >> ${OutfilePerBin}
		    #rm -fr _tmp_${refdir}_${dir}
		done # Close cycle over $bin
	        rm -fr _tmp_signal_all_${refdir}_${dir} _tmp_control_all_${refdir}_${dir} _tmp_all_${refdir}_${dir}
	    done # Close cycle over $state
	done # Close cycle over $binwidth
    done # Close cycle over $dir
done # Close cycle over $refdir
rm -fr _tmp_signal_all_${refdir}_${dir} _tmp_control_all_${refdir}_${dir} _tmp_all_${refdir}_${dir} _tmp_all_${refdir}_${dir}.R _tmp_${refdir}_${dir} _tmp_${refdir}_${dir}.R
sort -k 8d,8d -k 1d,1d -k 3d,3d -k 4n,4n ${OutfilePerBin} | uniq | awk '{if(NF==8) print $0}' > _a
mv _a ${OutfilePerBin}
awk '{printf("%-110s %s %-20s %s %s %s %s %s\n",$1,$2,$3,$4,$5,$6,$7,$8)}' ${OutfileOverall} > _a ; mv _a ${OutfileOverall}
