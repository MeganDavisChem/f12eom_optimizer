#!/bin/tcsh
#THE IS A SCRIPT FOR SUBMITTING F12+EOM JOBS ON ELAND COMPUTE NODES
#                        written by GST
#                        modified by TLE
#                        modified by MCD
#                   updated on Aug. 22, 2022


##
## SETUP ENVIRONMENTAL VARIABLES
##
#out of date
setenv TOP    /home/qc/
setenv QCTOP  $TOP/psi4conda
setenv QCPATH $QCTOP/bin
setenv SCRDRV /scratch

#activate conda environment with right version of psi4/optking/f12eom_optimizer
source /home/mdavis/anaconda3/etc/profile.d/conda.sh
conda activate p4env

setenv PATH "$PATH":/home/qc/molpro_2020.1/bin

#setenv PATH  "$QCPATH":/ddn/home1/r1621/maple/bin/tempQC/bin/cfour_v2devel_2015_06_18_g++_intel_2011.12.361/bin:/ddn/home4/r0683/QC/Bmrcc_2018-02-04:"$PATH"
#setenv LD_LIBRARY_PATH /ddn/home1/r1621/maple/bin/tempQC/bin/lib:/ddn/apps1/intel/lib/intel64
#setenv INTELPATH /usr/local/apps/intel/composer_xe_2011_sp1.12.361
#source $INTELPATH/mkl/bin/mklvars.csh intel64

#setenv PATH  "$QCPATH":/home/qc/cfour_v2devel_2015_06_18_g++_intel_2011.12.361/bin:"$PATH"
#setenv LD_LIBRARY_PATH /home/qc/psi4conda/ld_lib:
#rehash

setenv QC `which psi4`
if(! -e $QC) then
  echo "ERROR:  cannot find software installation"
  echo "        $QC"
  exit 1
endif

##
## SCRATCH SPACE
##
setenv WORKDIR `pwd`
if( $?PBS_JOBID ) then
  setenv PSI_SCRATCH $SCRDRV/$USER/PSI4_$PBS_JOBID
else
  set job_id = `date +%F_%T.%N`
  setenv PSI_SCRATCH $SCRDRV/$USER/PSI4_$job_id
endif
mkdir -p $PSI_SCRATCH
if($?) then
  echo "ERROR: unable to make $PSI_SCRATCH"
  exit
endif

##
## PRINT SOME USEFUL INFO 
##
echo " "
echo "Running threaded PSI4 on $HOST "
echo "$QC"
echo " "
echo "Working and scratch directories"
echo "$WORKDIR"
echo "$PSI_SCRATCH"
echo " "
if( $?NUM ) then
  setenv OMP_NUM_THREADS $NUM
  echo "Setting OMP_NUM_THREADS to $NUM"
  echo '(presumably from PBS file: setenv NUM $NCPUS)'
  echo "You can overide with -n flag"
  echo " "
else
  echo "use -n:  to specify number of threads"
  echo "use -i:  to specify input  file     (default = input.dat)"
  echo "use -o:  to specify output file     (default = file.out )"
  echo " "
  echo "-h or --help for more options"
  echo " "
endif

echo "Running..."
#echo "psi4 $argv"
f12eom_optimizer > output.dat
#psi4 $argv
echo "...Done"
echo " "

##
## CLEAN UP
##
cd $WORKDIR
rm -Rf $PSI_SCRATCH
