function r = psi_corrlen(fbase_in, fname_out, nfiles)
%
%   psi_corrlen(fbase_in, fname_out, nfiles)
%
%   normalize and average several correlation function and
%   find correlation length (half-width at half-height)
%
%   Input:
%       fbase_in    name base for input files
%       fname_out   text output file
%       nfiles      array with file numbers to be averaged 

  f = 0;


  %-- normalize and average correlation functions --

  for n = nfiles

    fname_in = [fbase_in, '.',  num2str(n,'%04d')];
    d = importdata(fname_in);

    f = f + d(:,2)/d(1,2);

  end

  f = f/length(nfiles);


  %-- correlation lenght -- 

  i1 = max(find(f>0.5));
  f1 = f(i1);
  f2 = f(i1+1);
    
  r = i1 + (0.5-f1)/(f2-f1) - 1;


  %-- output -- 

  %fprintf(' r =  %12.4e\n', r);

  fid = fopen(fname_out,'wt');

  fprintf(fid, '%% created by \"psi_corrlen.m\" from \"%s.*\"\n', fbase_in);
  fprintf(fid, '%% correlation length R = %12.4e (gridpoints)\n', r);
  fprintf(fid, '%% 1.radius (in gridpoints)  2.correlation function\n\n');

  for n=1:length(f)

    fprintf(fid, '%6d  %16.8e\n', (n-1), f(n));

  end

  fclose(fid);


end

%---------------------------------------------------------
