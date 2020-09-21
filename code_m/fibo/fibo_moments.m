function fibo_moments(fbase, fname_out, fnums)
%
%  "fibo_timeavg" is the script to compute time averages for occupation numbers
%  and ingegrals of motion over specified time intervals. 
%  Data input is from files "*.ak", data output to text file.
%
%  Input parameters:
%   fbase      string base for input files
%   fname_out  string for the file
%   fnums      array of input files to process
%
 
   fname = [fbase, '.param'];
  
   load(fname);  % ('fnum', 'P', 'G', 'm', 'runtype', 'dt', 'isave', 'nsave'); 

   M = m; 

   A = zeros(nsave*length(fnums), 2*M+1);

    i=0;

    for fnum = fnums 

        fname    = [fbase, '.',  num2str(fnum, '%04d'), '.ak'];

        fid = fopen(fname, 'rb');

        a = fread(fid, nsave*(2*M+1), 'double');
        a = reshape(a, [nsave, 2*M+1]);

        fclose(fid);

        A(i*nsave+1: (i+1)*nsave, :) = a;
  
        i = i+1;

    end

    
    %-- basic averaging ---

    A = reshape(A, [nsave*length(fnums), 2*M+1]);
    ntot = length(A);

    t = A(:,1);
    ind = (1:M)*2;
    A = A(:, ind)  + 1i*A(:, ind+1);   % A(ntot, M);

    Fi = fibonacci(M);
   
    n = A.*conj(A);
    nk = sum(n, 1)/ntot;
    m2k = sum(n.^2, 1)/ntot  ./ nk.^2 / 2  ;
    m3k = sum(n.^3, 1)/ntot  ./ nk.^3 / 6  ;
    clear n

    aa = A(:, 1:end-1) .* A(:, 2:end);
    s1k = sum(aa, 1)/ntot;
    fk  = sum(aa.*conj(aa), 1)/ntot; 

    aa = A(:, 1:end-1) .* conj(A(:, 2:end));
    s2k = sum(aa, 1)/ntot;

    s1k = [s1k,0];
    s2k = [s2k,0];
    fk  = [fk,0];

    clear aa

    nnk  = nk(1:end-1) .* nk(2:end);  
    ss1k = s1k(1:end-1) .* conj(s1k(2:end));  
    ss2k = s2k(1:end-1) .* conj(s2k(2:end));

    nnk  = [nnk,0];
    ss1k = [ss1k,0];
    ss2k = [ss2k,0];
 
    gk = fk - nnk - 2*real(ss1k) - 2*real(ss2k);

   %-- print out averages ---

   fid = fopen(fname_out, "w");
   fprintf(fid, "# Moments  \n");
   fprintf(fid, "# computed by \'fibo_moments.m\' for \'%s\'. \n", fbase);
   fprintf(fid, "# 1.i  2.Fi   3.nk  4.m2k  5.m3k   6.gk  7.fk  8.nnk  9.2re(ss1k)  10.2re(ss2k) \n\n");
   for k=1:M

     fprintf(fid, "%4d  %8d     %16.8e   %16.8e  %16.8e     %16.8e %16.8e    %16.8e  %16.8e  %16.8e\n", ...
             k, Fi(k),   nk(k), m2k(k), m3k(k),    gk(k), fk(k),   nnk(k), 2*real(ss1k(k)),  2*real(ss2k(k)) );
    
   end

   fclose(fid);
    
return

end

%---------------------

