function fibo_spectra(fbase, fname_out, fnums)
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
  
   load(fname);  % ('fnum', 'P', 'G', 'm', 'v', 'runtype', 'dt', 'isave', 'nsave'); 

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
    A = A(:, ind)  + 1i*A(:, ind+1);

    Fi = fibonacci(M);
   
    n = A.*conj(A);
    nk = sum(n, 1)/ntot;
    clear n

    A = [zeros(ntot,2), A];
    q = conj(A(:, 3:end)) .* A(:,2:end-1) .* A(:,1:end-2);
    qk = sum(q, 1)/ntot;
    
    clear A
    clear q

    %-- spectra and fluxes --

    Fip2 = fibonacci(M+2);
    Fip1 = Fip2(1:M+1);
    Fi   = Fip2(1:M);
 
    Fip2 = Fip1(2:end);
    Fip1 = Fip1(2:end);

    iF = fliplr(Fi);
    iFm1 = [iF(2:end),0];
    iFm2 = [iF(3:end),0,0];

    i = 1:M;

    F1 = Fi .* nk;
    F2 = Fip1 .* nk;
    FN = iF .* nk .* (-1).^i;

    V = Fi.^v;

    vm2 = [0, 0, V(1:end-2) ];
    qkp2 = [qk(3:end), 0, 0];   

    J = imag(qk) .* vm2 ; 
    K = real(qkp2) .* V * 2 ;

    Jp1 = [J(2:end), 0];
    Jp2 = [J(3:end), 0, 0];

    PI =   (Fip1 .* Jp1 +  Fi  .* Jp2) * 2;
    IP =  (-iFm2 .* Jp1 + iFm1 .* Jp2 ) .* (-1).^i  * 2;

   %-- print out averages ---

   fid = fopen(fname_out, "w");
   fprintf(fid, "# Averages for flux \n");
   fprintf(fid, "# computed by \'fibo_fluxavg.m\' for \'%s\'. \n", fbase);
   fprintf(fid, "# 1.i  2.Fi   3.F1  4.F2  5.FN    6.K  7.J    8.PI  9.IP \n\n");
   for k=1:M

     fprintf(fid, "%4d  %8d    %16.8e   %16.8e  %16.8e    %16.8e %16.8e   %16.8e  %16.8e\n", ...
             k, Fi(k),  F1(k), F2(k), FN(k),   K(k),  J(k),    PI(k), IP(k) );
    
  end

  fclose(fid);
    
return

end

%---------------------

