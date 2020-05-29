function twomode_mi(fbase_in)


%--- data ---  

    fname = [fbase_in, '.param'];

    load(fname);  % ('ntot', 'nPhi', 'nN1', 'nN2', 'max_N1', 'max_N2', 'g1', 'g2', 'p1', 'p2' );

    fid = fopen([fbase_in, '.P123'], 'rb');      % P123(nN1, nN2, nPhi, 'uint32');

    P123 = fread(fid, nN1*nN2*nPhi, 'uint32');

    fclose(fid);

    P123 = reshape(P123, nN1, nN2, nPhi);

	 
%-- normalise distributions to one  --
 
    P12 = squeeze(sum(P123,3));
    P1  = squeeze(sum(P12,2));
    P2  = squeeze(sum(P12,1));

    dn1 = max_N1/nN1; 
    dn2 = max_N2/nN2;
    dphi = 2*pi/nPhi;

    P123 = P123/ntot/(dn1*dn2*dphi);
    P1   = P1/ntot/dn1;
    P2   = P2/ntot/dn2;

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

    I123 =  S1 + S2 - (S123 - log2(2*pi));


    printf('%s  %d   %6.4f  %6.4f  %6.4f  %6.4f\n', fbase_in,  ntot,  S123, S1, S2, I123);
 


end

%---------------------

