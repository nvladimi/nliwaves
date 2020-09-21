function fibo_basics(fbase, fnums, ev)
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

    %-- read data, no rescaling ---

    A = reshape(A, [nsave*length(fnums), 2*M+1]);
    ntot = length(A);

    t = A(:,1);
    ind = (1:M)*2;
    A = A(:, ind)  + 1i*A(:, ind+1);

    Fi = transpose(fibonacci(M));
   
    aa = A.*conj(A);

    A = [A, zeros(ntot,2)];

    aaa = A(:, 1:end-2) .* A(:,2:end-1) .* conj(A(:,3:end));
    aaa = aaa + conj(aaa); 

    F0 = sum(aaa,2); 
    F1 = aa * Fi;
    F2 = aa(:, 2:end) * Fi(1:end-1);

    F0avg = sum(F0)/ntot;
    F1avg = sum(F1)/ntot;
    F2avg = sum(F2)/ntot;
 
    if (ev > 0)
      ind=(1:ev:length(t));
      ti=t(ind); F0=F0(ind); F1=F1(ind); F2=F2(ind); aa=aa(ind,:);
      er0=F0/F0(1)-1; er1=F1/F1(1)-1; er2=F2/F2(1)-1;

      figure(1);
      plot(ti, F0, '-ok', ti, F1, 'xr', ti, F2, '-b');
      set(gca, "fontsize", 20);  grid("on");

      figure(2);
      plot(ti, er0, '-ok', ti, er1, 'xr', ti, er2, '-b');
      set(gca, "fontsize", 20);  grid("on");

      figure(3);
      %plot(ti, aa, '-k');
      plot(ti, aa, '-k', ti, aa(:,1), '+r', ti, aa(:,2), 'xb', ti, aa(:,3), '-ok');
      set(gca, "fontsize", 20);  grid("on");

   end

    return

    %-- print out basics ---


%    printf("%s | %4.2f  %4.2f  %5.2e  %5.2e  %5.3f  %5.1f | ", fbase,  g1, g2, p1, p2, dt, ntot/1e6);
 
%    printf("%8.3f  %8.2e  %9.2e  %6.4f  %6.4f | ", chi, T, dT/T,  n1/T,  n2/T);
%    printf("%5.1f  %5.1f  %5.1f  %5.1f  %5.1f  %5.1f\n",   pi/sqrt(n1), pi/sqrt(4*n2),  1/g1, 1/g2,  n1/p1, n2/p2  );
 
    
return

end

%---------------------

