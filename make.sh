#!/bin/bash

dir=slicer/Slicer3/Applications/CLI/GroupWiseRegistration/ 
if [ -z $dir ]; then
  mkdir $dir
fi

cp GroupWiseRegistration/*  $dir
make -f slicer/Slicer3-build/Applications/CLI/GroupWiseRegistration/Makefile
