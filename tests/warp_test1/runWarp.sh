#!/bin/csh

#../../slicer3-6/Slicer3-build/lib/Slicer3/Plugins/WarpVolume --inputVolume ../G119_128.nrrd --warp G119_128_4_deformation.nii.gz --resultsDirectory ./
../../slicer3-6/Slicer3-build/lib/Slicer3/Plugins/WarpVolume --inputVolume /projects/schiz/3Tdata/case01009/diff/01009-dwi-filt-Ed.nhdr --warp ../G119_128_4_deformation.nii.gz --resultsDirectory ./
