function twomode_basics(fbase, fnums, ev)
 
   fname = [fbase, '.param'];
  
   load(fname);  % ('fnum', 'Gamma', 'Rflux', 'dt', 'isave', 'nsave'); 
 
    A = zeros(nsave*length(fnums), 5);

    i=0;

    for fnum = fnums 

        fname    = [fbase, '.',  num2str(fnum, '%04d'), '.a1a2'];

        fid = fopen(fname, 'rb');

        a = fread(fid, nsave*5, 'double');
        a = reshape(a, [nsave, 5]);

        fclose(fid);

        A(i*nsave+1: (i+1)*nsave, :) = a;
  
        i = i+1;

    end

    %-- read data, no rescaling ---

    A = reshape(A, [nsave*length(fnums), 5]);
    ntot = length(A);

    t = A(:,1); 
    A = A(:, 2:5);

    b1 = A(:,1) + 1i * A(:,2);
    b2 = A(:,3) + 1i * A(:,4);

    N1 = b1.*conj(b1);
    N2 = b2.*conj(b2);

    n1 = sum(sum(N1))/ntot;
    n2 = sum(sum(N2))/ntot;

    ind=(1:ev:length(t));
    ti=t(ind); N1i=N1(ind); N2i=N2(ind);
    plot(ti, N1i, '-r', ti, N2i, '-b');
    set(gca, "fontsize", 20);  grid("on")



    %-- print out basics ---


    g1  = - Gamma(1);
    g2  = - Gamma(3);
    p1  =   Rflux(1);
    p2  =   Rflux(3);

    T = (p1 + 2*p2)/2/(g1+g2);
    dT = p1/g1 - 2*p2/g2;
    chi = (g1 + g2)^2 /2/T;
    beta = 1 + 2*p2/p1;
    alpha = dT/T *g1*g2/(g1+g2)*beta;


    printf("%s | %4.2f  %4.2f  %5.2e  %5.2e  %5.3f  | ", fbase,  g1, g2, p1, p2, dt);
 
    printf("%7.3f  %5.2e  %4.2e  %6.4f  %6.4f | ", chi, alpha,  T,  n1/T,  n2/T);
    printf("%5.1f  %5.1f  %5.1f  %5.1f  %5.1f  %5.1f\n",   pi/sqrt(n1), pi/sqrt(4*n2),  1/g1, 1/g2,  n1/p1, n2/p2  );
 
    
return

end

%---------------------

