#!/bin/bash

name=pip_sec2_truth_test9_12-12-23

#swif2 cancel -delete -workflow ${name}

swif2 create -workflow ${name}

#N_jobs=25      #Number of jobs
N_jobs=2
workdir=/work/clas12/dir
gcard_dir=$workdir/gcards
gcard=$gcard_dir/clas12-default_pip_sec2.gcard
out_path=/output/dir
output_dir=$out_path/test_dir
NEVENTS=10

mkdir -p $output_dir
mkdir -p $output_dir/gemc
mkdir -p $output_dir/decoded
mkdir -p $output_dir/cooked
mkdir -p $output_dir/log

# Create jobs for simulating events
for ((i = 1; i <= $N_jobs; i++))
do
	swif2 add-job -workflow ${name} \
    		-constraint general \
    		-account clas12 \
    		-partition priority \
    		-disk 1gb \
    		-ram 1gb \
    		-time 1h \
		"source /group/clas12/packages/setup.sh && module load clas12 && gemc $gcard -N=$NEVENTS -USE_GUI=0 -OUTPUT='hipo, $output_dir/gemc/out_$i.hipo'"
done

swif2 run -workflow ${name}
###-output ${output_file_name}_$i.txt ${out_path} \
