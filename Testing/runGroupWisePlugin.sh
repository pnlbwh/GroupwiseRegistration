#!/bin/bash

data=../Data/Input
IMAGE1=$data/FA1.nii.gz
IMAGE2=$data/FA2.nii.gz
../GroupwiseRegistration --volumeFileNames $IMAGE1,$IMAGE2 --resultsDirectory ./ --num_multires_levels 5 --num_outer_iterations 10 --num_inner_iterations 40 --initial_sigma_diff 10 --final_sigma_diff 2 --reg_weight 10.0
