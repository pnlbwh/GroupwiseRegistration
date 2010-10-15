#!/bin/bash

#/Users/reckbo/Documents/slicers/trunk/Slicer3-build/lib/Slicer3/Plugins/GroupWiseRegistration --resultsDirectory ./ --outputVolume_Warp masked_case00747.nii.gz --num_multires_levels 6 --num_outer_iterations 10 --num_inner_iterations 50 --initial_sigma_diff 10 --final_sigma_diff 2 --reg_weight 10 --useJacFlag 
#/Users/reckbo/Documents/slicers/trunk/Slicer3-build/lib/Slicer3/Plugins/GroupWiseRegistration --volumeFileNames /Users/reckbo/Documents/slicers/trunk/Slicer3/Applications/CLI/GroupWiseRegistration/masked_case00747.nii.gz,/Users/reckbo/Documents/slicers/trunk/Slicer3/Applications/CLI/GroupWiseRegistration/masked_caseG119.nii.gz --resultsDirectory /Users/reckbo/Documents/slicers/trunk/Slicer3-build/Applications/CLI/GroupWiseRegistration --outputVolume_Warp masked_case00747.nii.gz --output_image output.nii.gz --num_multires_levels 6 --num_outer_iterations 10 --num_inner_iterations 50 --initial_sigma_diff 10 --final_sigma_diff 2 --reg_weight 10 --useJacFlag 
#./slicer/Slicer3-build/lib/Slicer3/Plugins/GroupWiseRegistration --volumeFileNames 10_masked_cases_affined/masked_case00747.nii.gz,10_masked_cases_affined/masked_caseG119.nii.gz --resultsDirectory ./ --outputVolume_Warp masked_case00747-norm.nii.gz --output_image output.nii.gz --num_multires_levels 6 --num_outer_iterations 10 --num_inner_iterations 50 --initial_sigma_diff 10 --final_sigma_diff 2 --reg_weight 10 --useJacFlag 
#/Users/reckbo/Documents/slicers/trunk/Slicer3-build/lib/Slicer3/Plugins/GroupWiseRegistration --volumeFileNames /Users/reckbo/Documents/slicers/trunk/Slicer3/Applications/CLI/GroupWiseRegistration/masked_case00747.nii.gz,/Users/reckbo/Documents/slicers/trunk/Slicer3/Applications/CLI/GroupWiseRegistration/masked_caseG119.nii.gz --resultsDirectory /Users/reckbo/Documents/slicers/trunk/Slicer3-build/Applications/CLI/GroupWiseRegistration --outputVolume_Warp masked_case00747.nii.gz --output_image slicer:0x10e2a79e0#vtkMRMLScalarVolumeNode1 --num_multires_levels 6 --num_outer_iterations 10 --num_inner_iterations 50 --initial_sigma_diff 10 --final_sigma_diff 2 --reg_weight 10 --useJacFlag 

#./slicer3-6/Slicer3-build/lib/Slicer3/Plugins/GroupWiseRegistration --volumeFileNames 10_masked_cases_affined/masked_case00747.nii.gz,10_masked_cases_affined/masked_caseG119.nii.gz --resultsDirectory ./ --outputVolume_Warp masked_case00747-norm.nii.gz --output_image output.nii.gz --num_multires_levels 5 --num_outer_iterations 1 --num_inner_iterations 5 --initial_sigma_diff 10 --final_sigma_diff 2 --reg_weight 10 --useJacFlag 
#./slicer3-6/Slicer3-build/lib/Slicer3/Plugins/GroupWiseRegistration --volumeFileNames 10_masked_cases_affined/masked_caseG139.nii.gz,10_masked_cases_affined/masked_caseG137.nii.gz,10_masked_cases_affined/masked_case0581.nii.gz,10_masked_cases_affined/masked_case0628.nii.gz,10_masked_cases_affined/masked_case0652.nii.gz,10_masked_cases_affined/masked_case00747.nii.gz,10_masked_cases_affined/masked_caseG119.nii.gz --resultsDirectory ./ --outputVolume_Warp masked_case00747-norm.nii.gz --output_image output.nii.gz --num_multires_levels 5 --num_outer_iterations 1 --num_inner_iterations 5 --initial_sigma_diff 10 --final_sigma_diff 2 --reg_weight 10 --useJacFlag 
#./slicer3-6/Slicer3-build/lib/Slicer3/Plugins/GroupWiseRegistration --volumeFileNames 10_masked_cases_affined/masked_case0628.nii.gz,10_masked_cases_affined/masked_case0652.nii.gz,10_masked_cases_affined/masked_case00747.nii.gz,10_masked_cases_affined/masked_caseG119.nii.gz --resultsDirectory ./ --outputVolume_Warp masked_case00747-norm.nii.gz --output_image output.nii.gz --num_multires_levels 5 --num_outer_iterations 1 --num_inner_iterations 5 --initial_sigma_diff 10 --final_sigma_diff 2 --reg_weight 10 --useJacFlag 
#./slicer3-6/Slicer3-build/lib/Slicer3/Plugins/GroupWiseRegistration --volumeFileNames 10_masked_cases_affined/masked_case00747.nii.gz,10_masked_cases_affined/masked_caseG119.nii.gz --resultsDirectory ./ --outputVolume_Warp masked_case00747-norm.nii.gz --output_image output.nii.gz --num_multires_levels 5 --num_outer_iterations 1 --num_inner_iterations 5 --initial_sigma_diff 10 --final_sigma_diff 2 --reg_weight 10 --useJacFlag 
#./slicer3-6/Slicer3-build/lib/Slicer3/Plugins/GroupWiseRegistration --volumeFileNames 10_masked_cases_affined/masked_case0628.nii.gz,10_masked_cases_affined/masked_case0652.nii.gz,10_masked_cases_affined/masked_caseG119.nii.gz --resultsDirectory ./ --outputVolume_Warp masked_case00747-norm.nii.gz --output_image output.nii.gz --num_multires_levels 5 --num_outer_iterations 2 --num_inner_iterations 5 --initial_sigma_diff 10 --final_sigma_diff 2 --reg_weight 10 --useJacFlag 
#./slicer3-6/Slicer3-build/lib/Slicer3/Plugins/GroupWiseRegistration --volumeFileNames 10_masked_cases_affined/masked_case0652.nii.gz,10_masked_cases_affined/masked_caseG119.nii.gz --resultsDirectory ./ --outputVolume_Warp masked_case00747-norm.nii.gz --output_image output.nii.gz --num_multires_levels 5 --num_outer_iterations 2 --num_inner_iterations 2 --initial_sigma_diff 10 --final_sigma_diff 2 --reg_weight 10 --useJacFlag 
#./slicer3-6/Slicer3-build/lib/Slicer3/Plugins/GroupWiseRegistration --volumeFileNames ./input/00747_128.nii.gz,./input/0652_128.nii.gz,./input/g119_128.nii.gz --resultsDirectory ./output --output_image output.nii.gz --num_multires_levels 5 --num_outer_iterations 5 --num_inner_iterations 5 --initial_sigma_diff 10 --final_sigma_diff 2 --reg_weight 10 --useJacFlag 
#./slicer3-6/Slicer3-build/lib/Slicer3/Plugins/GroupWiseRegistration --volumeFileNames ./input/0652_128.nii.gz,./input/g119_128.nii.gz --resultsDirectory ./output --output_image output.nii.gz --num_multires_levels 5 --num_outer_iterations 5 --num_inner_iterations 5 --initial_sigma_diff 10 --final_sigma_diff 2 --reg_weight 10 --useJacFlag 


./slicer3-6/Slicer3-build/lib/Slicer3/Plugins/GroupWiseRegistration --volumeFileNames ./input_128/0652_128.nii.gz,./input_128/g119_128.nii.gz --resultsDirectory ./output --output_image output.nii.gz --num_multires_levels 5 --num_outer_iterations 5 --num_inner_iterations 5 --initial_sigma_diff 10 --final_sigma_diff 2 --reg_weight 10 --useJacFlag 


