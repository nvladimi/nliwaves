function twomode_array_n1n2theta(fbase, fnums, fbaseout)

tic
    nPhi    = 32;
    nN1     = 200;
    nN2     = 100;
    max_N1  = 20;
    max_N2  = 10;

    fname = [fbase, '.param'];
    load(fname);  % ('fnum', 'Gamma', 'Rflux', 'dt', 'isave', 'nsave'); 
      
    ntot = nsave*length(fnums)

    B = zeros(ntot, 5);

    i=0;

    for fnum = fnums 

        fname    = [fbase, '.',  num2str(fnum, '%04d'), '.a1a2'];

        fid = fopen(fname, 'rb');

        b = fread(fid, nsave*5, 'double');
        b = reshape(b, [nsave, 5]);

        fclose(fid);

        B(i*nsave+1: (i+1)*nsave, :) = b;
  
        i = i+1;

    end

    
   % B = B(1:ntmp,:);
   % ntot = ntmp;    

    %-- rescale data -- 

    g1  = - Gamma(1);
    g2  = - Gamma(3);
    p1  =   Rflux(1);
    p2  =   Rflux(3);

    T = (p1 + 2*p2)/2/(g1+g2);
    dT = p1/g1 - 2*p2/g2;
    chi = (g1 + g2)^2 /2/T;
    beta = 1 + 2*p2/p1;
    alpha = dT/T *g1*g2/(g1+g2)*beta;

    g = 0.5*g1*g2/(g1+g2);

    B = reshape(B, [ntot, 5]);
    t = B(:,1) * g ;
    B = B(:, 2:5) /sqrt(T);

    %-- convert data to probability variables -- 

    b1 = B(:,1) + 1i * B(:,2);
    b2 = B(:,3) + 1i * B(:,4);

    %F =  4 * imag ( b1 .* b1 .* conj(b2) );
    %K =  2 * real ( b1 .* b1 .* conj(b2) );

    N1 = b1.*conj(b1);
    N2 = b2.*conj(b2);

    N1_avg = sum(N1)/ntot
    N2_avg = sum(N2)/ntot

    a1=angle(b1);
    a2=angle(b2);
    phi = 2*a1 - a2;
    ind = find(phi > pi);  phi(ind) = phi(ind) - 2*pi;
    ind = find(phi < -pi); phi(ind) = phi(ind) + 2*pi;

    %-- create probabilty array ---

  
    P123 = zeros(nN1, nN2, nPhi+1, 'uint32');
    nPi = nPhi/2;


    %-- convert to index and add counts---
    
    iN1 = floor( N1/max_N1 * nN1 ) + 1;
    iN2 = floor( N2/max_N2 * nN2 ) + 1;
    iN1  = min(iN1, nN1);
    iN2  = min(iN2, nN2);

    iPhi = round( phi./pi * nPi) + nPi + 1;
toc

tic

  for i=1:ntot
 	    P123(iN1(i),iN2(i),iPhi(i)) = P123(iN1(i),iN2(i),iPhi(i))  + uint32(1);
   end

toc

   P123(:,:,1) = P123(:,:,1) + P123(:,:,end);
   P123 = P123(:,:,1:end-1);
   
   fnameout = [fbaseout, '.param'];


save( fnameout, 'ntot', 'nPhi', 'nN1', 'nN2', 'max_N1', 'max_N2', 'g1', 'g2', 'p1', 'p2' ); 

fid = fopen([fbaseout, '.P123'], 'wb');
fwrite(fid, P123, 'uint32');
fclose(fid);

return


end

%---------------------
