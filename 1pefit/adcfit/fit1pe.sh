#!/bin/bash
## Example
##   make
##   ./fit1pe.sh inputsample.root

inputfile=$1
#npeak=2
npeak=3

./performance1pe $inputfile ${inputfile}_output.root $npeak 0 

