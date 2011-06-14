#!/bin/csh

../../slicer3-6/Slicer3-build/lib/Slicer3/Plugins/GroupWiseRegistration --volumeFileNames ../G119_128.nrrd,../G132_128.nrrd,../G133_128.nrrd,../G136_128.nrrd,../G137_128.nrrd,../G139_128.nrrd --resultsDirectory ./ --output_image atlas.nii.gz --num_multires_levels 5 --num_outer_iterations 5 --num_inner_iterations 5 --initial_sigma_diff 10 --final_sigma_diff 2 --reg_weight 10.0
