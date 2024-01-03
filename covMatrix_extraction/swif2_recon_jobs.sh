#!/bin/bash

name=pip_recon_test3_12-19-23

#swif2 cancel -delete -workflow ${name}

swif2 create -workflow ${name}

#workdir=/work/clas12/dir
#out_path=/output/dir
#output_dir=$out_path/test_dir

#input_dir=/path/to/gemc/directory
#input_dir=/volatile/clas12/reedtg/clas12_kinfitter/cov_matrix/pip_sec2_truth/gemc
output_dir=/volatile/clas12/reedtg/clas12_kinfitter/cov_matrix/pip_sec2_test
input_dir=$output_dir/gemc


file_count=0
max_files=2
# Create jobs for reconstructing simulation
for input_file in $input_dir/*
do
    # Uncomment to break loop after some number of files are read
    if [ $file_count -ge $max_files ]; then
        break  
    fi

    # Extracting the filename without extension
    filename=$(basename -- "$input_file")
    base_filename="${filename%.*}"
	swif2 add-job -workflow ${name} \
    		-constraint general \
    		-account clas12 \
    		-partition priority \
    		-disk 1gb \
    		-ram 2gb \
    		-time 1h \
		"source /group/clas12/packages/setup.sh && module load clas12 && \ 
			$COATJAVA/bin/recon-util -i $input_dir/$base_filename.hipo -o $output_dir/cooked/$base_filename.rec.hipo -c 2  > $output_dir/log/rec_$base_filename.log"
    ((file_count++))
done

swif2 run -workflow ${name}
###-output ${output_file_name}_$i.txt ${out_path} \
