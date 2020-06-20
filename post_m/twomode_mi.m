function twomode_mi(fbase)
%
%  reads parameter file for probability distribution (*_post??.param)
%  reads propability distribution P(n1,n2,theta)
%  computes MI and entropies for this distributon
% 
%  header
%  #1.fname  2.samples(M)  3.chi  4.dT/T   5.S12    6.S1     7.S2     8.S3    9.I12
%

%--- data ---  

    fname = [fbase, '.param'];

    load(fname);  % ('ntot', 'nPhi', 'nN1', 'nN2', 'max_N1', 'max_N2', 'g1', 'g2', 'p1', 'p2' );

    fid = fopen([fbase, '.P123'], 'rb');      % P123(nN1, nN2, nPhi, 'uint32');

    P123 = fread(fid, nN1*nN2*nPhi, 'uint32');

    fclose(fid);

    P123 = reshape(P123, nN1, nN2, nPhi);

    P1  = squeeze(sum( sum(P123,3), 2));
    P2  = squeeze(sum( sum(P123,3), 1));
    P3  = squeeze(sum( sum(P123,2), 1));
	 
%-- normalise distributions to one  --

    dn1 = max_N1/nN1; 
    dn2 = max_N2/nN2;
    dphi = 2*pi/nPhi;

    P123 = P123/ntot/(dn1*dn2*dphi);
    P1   = P1/ntot/dn1;
    P2   = P2/ntot/dn2;
    P3   = P3/ntot/dphi;

    %figure(1)
    %imagesc(log(rot90(P12))); axis("equal", "off"); colormap("jet");

    %figure(2)
    %semilogy(1:length(P1), P1, '-or', 1:length(P2), P2, '-ob'); 

    
%-- compute mutual information  --

    q =  reshape(P123, nN1*nN2*nPhi, 1);

    ind=find(q>0);
    S123 = - sum(q(ind).*log2(q(ind))) * (dn1*dn2*dphi);
  
    ind=find(P1>0);
    S1 = - sum(P1(ind).*log2(P1(ind))) * dn1;

    ind=find(P2>0);
    S2 = - sum(P2(ind).*log2(P2(ind))) * dn2;

    ind=find(P3>0);
    S3 = - sum(P3(ind).*log2(P3(ind))) * dphi;


    I12 =  S1 + S2 + log2(2*pi) - S123;


    warning("off");
    T = (p1 + 2*p2)/2/(g1+g2);
    dT = p1/g1 - 2*p2/g2;
    chi = (g1 + g2)^2 /2/T;

 
    printf('%s  %4d  %8.3f  %7.4f  %7.4f  %7.4f  %7.4f  %7.4f  %7.4f\n', ...
	   fbase,  ntot/1e6,  chi,  dT/T,  S123, S1, S2, S3, I12);
 
end

%---------------------

