#!/bin/sh
#
# Script for submitting up to 8 serial jobs on a single node of
# Pequena.  The script is meant to be called from a PBS submission 
# script, for example "pequena_serial.sub".  Executable "getid.x" 
# has to be compiled separately.
#
# Important: change permissions on this file to 'rwx'.

. ${MODULESHOME}/init/sh
module load  fftw/3.2.1/intel

myid=`$HOME/work/bin/getid.x`
program=$HOME/work/bin/kelseg.x

cd $HOME/scratch/kelseg

case $myid
in
  0) $program -f m41a -A 1.025 -N0 256 -L 25.6 -Ngrids 6 ;;
  1) $program -f m41b -A 1.025 -N0 512 -L 51.2 -Ngrids 6 ;;
  2) $program -f m42a -A 1.050 -N0 256 -L 25.6 -Ngrids 6 ;;
  3) $program -f m42b -A 1.050 -N0 512 -L 51.2 -Ngrids 6 ;;
esac
