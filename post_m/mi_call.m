#!/opt/octave/bin/octave -qf

fbase ='../test/a5';

fbaseout = '../test/test02';

m1 = 6;
m2 = 12;

ndiag  = 1000;               % number of samples per file;
fnums  = [30:31];            % array of files numbers to process

amp0      = 0.15;            % max amplitude of a mode
nAmp      =   64;            % number of bins for amplitudes
nPhi      =   16;            % number of bins for phase (must be even!)

N         =   80;            % saved number of modes 
N0        =  256;            % total number of modes

%--------------------------------

mi_core(fbase, fbaseout, fnums, N, N0, ndiag, amp0, nAmp, nPhi, m1, m2);

%--------------------------------
