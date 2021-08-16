function fibo_probnk(fbasein, fbaseout, fstart)
%
%  "fibo_basics" is the script to compute derived parameteres of the run
%  and occupation numbers, and to estimate time step restrictions. 
%  Data input is from files "*.ak", data output is to the screen.
%
%  Input parameters:
%   fbase     string base for input files
%   fnums     array of input files to process
%   ev        if ev>0, plot occupation numbers using every "ev" datapoint
%

   nbins = 200;       % number of bins for each mode
   tail_factor = 40;  % extend of bins above average

   seeds = [1:4];      % realizations 

%-- first round: read data and compute average --

   ntot = 0;

   for s = seeds 

       fbase = [fbasein, '_s', num2str(s)];
       fbase = [fbase, '/', fbase];

       load([fbase, ".param"]);

       M = m; 

       if (s == 1) 
          aa_avg = zeros(1,M);
       end

       fnum = fstart;

       while 1

           fname    = [fbase, '.',  num2str(fnum, '%04d'), '.ak'];

           if exist(fname, "file")

              fid = fopen(fname, 'rb');
              a = fread(fid, nsave*(2*M+1), 'double');
              a = reshape(a, [nsave, 2*M+1]);

              ind = (1:M)*2;
              a = a(:, ind)  + 1i*a(:, ind+1);

              aa = a.*conj(a);  % <-- aa(nsave,M)

              aa_avg = aa_avg + sum(aa,1);
              ntot   = ntot + nsave;
              fnum   = fnum + 1;

           else % file does not exist
             break
           end

       end  % loop over files

   end % loop over seeds

   aa_avg = aa_avg/ntot;

%-- debugging --

   fnameout = [fbaseout, '.param'];

   save( fnameout, 'fbase', 'M', 'ntot', 'nbins', 'tail_factor', 'aa_avg' );

%-- create probability array  --

    Paa = zeros(nbins, M, 'uint32');

%-- second round: read data and assemble distribution --

   for s = seeds

       fbase = [fbasein, '_s', num2str(s)];
       fbase = [fbase, '/', fbase];
 
       load([fbase, ".param"]);

       fnum = fstart;

       while 1
 
           fname    = [fbase, '.',  num2str(fnum, '%04d'), '.ak'];

           if exist(fname, "file")

              fid = fopen(fname, 'rb');
              a = fread(fid, nsave*(2*M+1), 'double');
              a = reshape(a, [nsave, 2*M+1]);

              ind = (1:M)*2;
              a = a(:, ind)  + 1i*a(:, ind+1);

              aa = a.*conj(a);  % <-- aa(nsave,M)

              for m=1:M

                 max_bin = tail_factor * aa_avg(m);
                 aam = aa(:,m);
                 iaa = floor( aam/max_bin * nbins ) + 1;
                 iaa  = min(iaa, nbins);

                 for i=1:nsave
                      Paa(iaa(i),m) = Paa(iaa(i),m)  + uint32(1);
                 end
              end % loop over modes

              fnum = fnum+1;
 
           else % file does not exist
             break
           end

       end  % loop over files

   end % loop over seeds


%---------------------------------------
   
   fnameout = [fbaseout, '.param'];

   save( fnameout, 'fbase', 'M', 'ntot', 'nbins', 'tail_factor', 'aa_avg' ); 

   fid = fopen([fbaseout, '.Paa'], 'wb');
   fwrite(fid, Paa, 'uint32');
   fclose(fid);

return

end

%---------------------

