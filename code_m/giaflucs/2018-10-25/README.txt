--------------------------
Files in current directory
--------------------------

giaflucs.m                 core computing, propagation
giaflucs_profiles.m        run/rerun selected seeds, compute profiles and FWHM
giaflucs_screen.m          generates single screen (former one_screen_smooth.m)
run_mseed.m                supplementary script used in multiseed runs         
run_mset.m                 supplementary script used in multiset runs

legacy.param               example of input parameters
multiseed.param            example of input parameters for multiseed runs
set001.param               example of input parameters for multiset runs 
set002.param               example of input parameters for multiset runs
set003.param               example of input parameters for multiset runs

set001_seeds.dat           example file with seed for profile averaging

comet_multiseed.sub        submission script specific to comet
comet_multiset.sub         submission script specific to comet
stampede_multiseed.sub     submission script specific to stampede2
stampede_multiset.sub      submission script specific to stampede2
odin_multiseed.sh          submission script specific to odin
odin_multiset.sh           submission script specific to odin

README.txt                this file

-----------------------------
Code modification and running
-----------------------------

To run the code, modify only *.param, *.sub, and *.sh files.
All matlab/octave files are supposed to be used without modification.
Be cautious applying old scripts to new data, as data format
might have changed.

The dataflow structure is the same for all three machines. 
Job submission:

stampede2>   sbatch stampede_multiseed.sub
comet>       sbatch    comet_multiseed.sub
odin>                 ./odin_multiseed.sh
or
stampede2>   sbatch  stampede_multiset.sub
comet>       sbatch     comet_multiset.sub
odin>                  ./odin_multiset.sh

To run on a personal computer:

octave> giaflucs(fbase, seed);
or 
octave> giaflucs_profiles(fbase);

The code will read parameters from corresponding "fbase.param" file
which must to be prepared before job submission.

Login to Comet as username@comet.sdsc.xsede.org with XSEDE password.
Login to Stampede2 as username@stampede2.tacc.utexas.edu with TACC 
password and token.

It is more efficeint to run "multiseed" setup on Stampede2 and 
"multiset" setup on Comet.

----------------------------------------
Types of simulations: (i) multiseed mode
----------------------------------------

This mode is intended for comuting large ensemble of seeds without 
profile averaging.  The tasks are parallelized over seeds - that is
the first task computes seeds [1:n], the second task computes seeds 
[n+1:2n] and so on.  All simulations use the same parameter file
"fbase.param".  See "multiseed.param" as an example. The basename of
the files can be arbitrary. The second part of the parameter file,
which specifies parameters for profile averaging, is ignored.
The code creates multiple files with maximum intensity versus
propagarion distance, "fbase.00001', etc. They are placed in the 
directory specified in the parameter file; make sure this directory
exists. Since parallelization over seeds is straightforward to
implement, multiseed mode is more efficient to run on Stampede2.

----------------------------------------
Types of simulations: (ii) multiset mode
----------------------------------------

This mode is intended for computing multiple sets with different 
simulation parameters on the same computing node. The code uses 
a different parameter file for each set. They must be named in 
the particular form,  "set001.param", "set002.param", etc.

Averaged profile is computed for each set. Parameters for averaging
are specified in the second half of parameter file. Notice that
we have an option to run seeds preselected in "set001_seeds.dat" etc,
or to run all seeds consecutively.  In the second case, averaging 
can be done for selected interval of intensity (with byproduct 
of "set001_seeds.dat", etc). The code writes the averaged profiles, 
I(r) with FWHM in the header to  to "set001_ravg.dat", etc, in addition 
to "fbase.00001', etc, in specified directory (make sure it exists).

In multiset mode, faster cores are preferred to large number of 
cores. It is unreasonable to write 272 parameters files to make 
small datasets on Stampede, as opposed to 24 parameter files for 
larger datasets on Comet.  

In the test setups provided, "set001" is computing averaged profile
for 50 preselected seed, "set002" is running 200 seeds with 
10 screens and "set003" is running 100 seeds with 20 screens. 
In both "set002" and "set003" about 10 seeds are selected for 
profile averaging. This job takes 18 min on Comet, and we can add 
22 more setups to it and extend it for larger ensembles.

----------------------------------------


