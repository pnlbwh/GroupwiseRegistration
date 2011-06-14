#!/bin/bash

DATA_DIR='/net/lmi-terastation-b/mnt/array1/share/pnl/home/sylvain/data/ControlsAtlas/FAinAtlas'
IMAGE1=$DATA_DIR/01183-Rgd-fa-Atlas.nii.gz
IMAGE2=$DATA_DIR/01181-Rgd-fa-Atlas.nii.gz
#../src/GroupWiseRegistration --volumeFileNames $DATA_DIR/G119_128.nrrd,$DATA_DIR/G132_128.nrrd --resultsDirectory ./ --num_multires_levels 5 --num_outer_iterations 5 --num_inner_iterations 5 --initial_sigma_diff 10 --final_sigma_diff 2 --reg_weight 10.0
../src/GroupWiseRegistration --volumeFileNames $IMAGE1,$IMAGE2 --resultsDirectory ./ --num_multires_levels 5 --num_outer_iterations 5 --num_inner_iterations 5 --initial_sigma_diff 10 --final_sigma_diff 2 --reg_weight 10.0
