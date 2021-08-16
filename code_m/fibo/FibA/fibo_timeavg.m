function fibo_timeavg(fbase, fname_out, setsize, fnums)
%
%  "fibo_timeavg" is the script to compute time averages for occupation numbers
%  and ingegrals of motion over specified time intervals. 
%  Data input is from files "*.ak", data output to text file.
%
%  Input parameters:
%   fbase      string base for input files
%   fname_out  string for the file
%   setsize    length of time interval, in terms of number of samples
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

    
    %-- read data, no rescaling ---

    A = reshape(A, [nsave*length(fnums), 2*M+1]);
    ntot = length(A)
    numsets = ntot/setsize

    t = A(:,1);
    ind = (1:M)*2;
    A = A(:, ind)  + 1i*A(:, ind+1);

    Fi = transpose( fibonacci(M) );
    iF = transpose( fliplr(fibonacci(M)) .* (-1).^(1:M) );

    V = Fi.^v;
   
    aa = A.*conj(A);
    
    A = [A, zeros(ntot,2)];

    aaa = A(:, 1:end-2) .* A(:,2:end-1) .* conj(A(:,3:end));
    aaa = (aaa + conj(aaa)); 

    F0 = aaa * V; 
    F1 = aa * Fi;
    F2 = aa(:, 2:end) * Fi(1:end-1);
    FN = aa * iF;

    %-- averaging --

    aa_avg  = zeros(numsets, M);
    F0_avg = zeros(numsets, 1);
    F1_avg = zeros(numsets, 1);
    F2_avg = zeros(numsets, 1);
    FN_avg = zeros(numsets, 1);
    t_avg = zeros(numsets, 1);

   for k=1:numsets
	ind = (1:setsize) + setsize*(k-1);
        aa_avg(k,:) = sum(aa(ind,:),1)/setsize;
        F0_avg(k) =  sum(F0(ind))/setsize;
        F1_avg(k) =  sum(F1(ind))/setsize;
        F2_avg(k) =  sum(F2(ind))/setsize;
        FN_avg(k) =  sum(FN(ind))/setsize;
        t_avg(k)  =  sum(t(ind))/setsize;
   end

   out_avg = [t_avg, aa_avg, F0_avg, F1_avg, F2_avg,  FN_avg];


   %-- print out averages ---

   fid = fopen(fname_out, "w");
   fprintf(fid, "# Averages for occupation number and integerals \n");
   fprintf(fid, "# computed by \'fibo_timeavg.m\' for \'%s\' using setsize %d. \n", fbase, setsize);
   fprintf(fid, "# 1.t  2.-%d.|a_k|^2  %d.<aaa*+cc>  %d.F1  %d.F2  %d.FN\n\n", M+1,  M+2, M+3, M+4, M+5);

   for k=1:numsets

     fprintf(fid, "%16.8e  ", out_avg(k,:));
     fprintf(fid, "\n");
    
  end

  fclose(fid);
    
return

end

%---------------------

