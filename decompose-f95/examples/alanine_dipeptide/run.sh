#!/bin/tcsh

# clean out output directory
rm -rf output
mkdir output

# run state decomposition
../../src/decompose control-test.xml

