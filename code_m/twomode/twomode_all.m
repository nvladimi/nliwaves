function twomode_all()
%
% cycles for files in a dataset to compute and print to screen MI and entropies; 
% calls "twomode_mi" by default; can work with "twomode_mi_theory" or "twomode_basics" 
%

%-- header for twomode_mi(fname)
%   #1.fname  2.samples(M)  3.chi  4.dT/T   5.S12    6.S1     7.S2     8.S3    9.I12


  prefix='POST/Post02/' ;
  suffix='_post02' ;
  name1 = ['x1'; 'x23'; 'x27'; 'x3'; 'x32'; 'x36'; 'x4'; 'x40'; 'x43'; 'x5'; 'x6'; 'x7'];
  name2 = ['d1'; 'd2'; 'd3'; 'd4'; 'd5'; 'd6';  'd7'];

  %prefix='POST/Post06/' ;
  %suffix='_post06' ;
  %name1 = ['x2'; 'x23'; 'x27'; 'x3'; 'x32'; 'x36'; 'x40'; 'x43'; 'x5'; 'x6'];
  %name2 = ['d1'; 'd2'; 'd3'; 'd4'; 'd5'; 'd6';  'd7'];

  %prefix='POST/Post03/' ;
  %suffix='_post03' ;
  %name1 = ['y02'; 'y03'; 'y04'; 'y05'; 'y08'; 'y10'; 'y13'; 'y20'];
  %name2 = ['pA'; 'p0'; 'p1'; 'p2'; 'p3'; 'p4'; 'p5'; 'p6';  'p7';  'p8';  'p9'; 'pB'];



%-- core cycle, use with on the the following: --
%-- twomode_mi(fname),  twomode_mi_theory (fname),  twomode_basics(fname, fnums, ev) --

  printf("\n");

  for i1=1:length(name1)
     for i2=1:length(name2)
	  fname = [prefix, deblank(name1(i1,:)), deblank(name2(i2,:)), suffix];
          twomode_mi(fname);
      end
      printf("\n");
   end;

 
end

%-----------------------------
