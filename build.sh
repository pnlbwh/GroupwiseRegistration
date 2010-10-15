#!/bin/bash

dir="../slicer3-6/Slicer3/Applications/CLI/GroupWiseRegistration"

if [ ! -e $dir ]; then
  mkdir $dir
  echo copying to $dir
  cp GroupWiseRegistration/*  $dir
  echo 'add_subdirectory(GroupWiseRegistration)' >> $dir/../CMakeLists.txt
  slicer3-6/Slicer3/Scripts/getbuildtest.tcl -t "" --no-slicer-update
else
  echo copying to $dir
  cp GroupWiseRegistration/*  $dir
fi

make -f slicer3-6/Slicer3-build/Applications/CLI/GroupWiseRegistration/Makefile
