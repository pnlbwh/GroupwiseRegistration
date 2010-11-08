#!/bin/bash

dir="slicer3-6/Slicer3/Applications/CLI/GroupWiseRegistration"

echo copying to $dir
cp -R WarpVolume $dir
cp CMakeLists.txt $dir

make -f slicer3-6/Slicer3-build/Applications/CLI/GroupWiseRegistration/Makefile
