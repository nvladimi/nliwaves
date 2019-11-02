function reavg_any(fbase_in, fname_out, nfiles)
%
%   reavg_any(fbase_in, fname_out, nfiles)
%
%   average several text files of the same format, header omitted
%
%   Input:
%       fbase_in    name base for input files
%       fname_out   text output file
%       nfiles      array with file numbers to be averaged 

  data = 0;

  for n = nfiles

    fname_in = [fbase_in, '.',  num2str(n,'%04d')];

    d = importdata(fname_in);
    % d(:,2)=d(:,2)/d(1,2);      %%% uncomment to normalize data (HACK!)
    data = data + d;

  end

  data = data/length(nfiles);

  [nmax, nfields] = size(data);


  %-- output -- 

  fid = fopen(fname_out,'wt');

  fprintf(fid, '%% created by \"reavg_all.m\" from \"%s.*\"\n', fbase_in);
  fprintf(fid, '%% header omitted\n\n');


  for n=1:nmax

    fprintf(fid,' %16.8e', data(n,:));
    fprintf(fid, '\n');

  end

  fclose(fid);


end

%---------------------------------------------------------
