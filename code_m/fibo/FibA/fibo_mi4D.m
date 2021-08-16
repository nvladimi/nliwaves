function fibo_mi4D(fbase, saveimages)
%
%  reads parameter file for probability distribution
%  reads propability distribution P4D(n1,n2,n3,theta)
%  computes MI and entropies for this distributon
% 
%  header
%  
%  #1.i0   2.MI      3.S4D       4.S1     5.S2     6.S3     7.S4   8.fname  9.samples(M) 
%

%--- data ---  

    fname = [fbase, '.param'];
    load(fname);  % ('i0', 'nBins', 'dB', 'ntot',  'aa_avg'); 

    fname = [fbase, '.P4D'];
    fid = fopen(fname, 'rb');      % P4(nN1, nN2, nN3, nPhi, 'uint32');

    P4D = fread(fid, nBins(1)*nBins(2)*nBins(3)*nBins(4), 'uint32');

    fclose(fid);
 
    P4D = reshape(P4D, nBins(1), nBins(2), nBins(3), nBins(4));


    P1  = squeeze(sum( sum( sum(P4D, 4), 3), 2));
    P2  = squeeze(sum( sum( sum(P4D, 4), 3), 1));
    P3  = squeeze(sum( sum( sum(P4D, 4), 2), 1));
    P4  = squeeze(sum( sum( sum(P4D, 3), 2), 1));

    P12  = squeeze(sum( sum(P4D,4), 3));
    P13  = squeeze(sum( sum(P4D,4), 2));
    P23  = squeeze(sum( sum(P4D,4), 1));

    % disp([sum(P1), ntot]);
	 
%-- normalise distributions to one  --

    dBins = aa_avg*dB;

    dphi = 2*pi/nBins(4);

    P4D  = P4D/ntot/(dBins(1)*dBins(2)*dBins(3)*dphi);
    P1   = P1/ntot/dBins(1);
    P2   = P2/ntot/dBins(2);
    P3   = P3/ntot/dBins(3);
    P4   = P4/ntot/dphi;

    P12  = P12/ntot/dBins(1)/dBins(2);
    P13  = P13/ntot/dBins(1)/dBins(3);
    P23  = P23/ntot/dBins(2)/dBins(3);

if saveimages

    warning('off');

    %l0=-5; levs =(l0:1:l0+15);
    %l0=-15; levs =(l0:1:10+15);
    l0=-12; levs =(l0:1:5);


    % figure(1);  imagesc(log(rot90(P13))); axis("square", "off"); colormap("jet");
    % figure(2);  imagesc(log(rot90(P23))); axis("square", "off"); colormap("jet");
    % figure(3);  imagesc(log(rot90(P12))); axis("square", "off"); colormap("jet");


figure(1); clf; q = max(log(P13), l0);
    contourf(flipud(rot90(q)),levs); axis("square"); colormap("jet");
    xlabel("0 < n( i-2 ) < 32 <n( i-2 )>");  ylabel("0 < n(i) < 32 <n(i)>");
    print([fbase,'_lnP13_lev.png'], '-dpng');

figure(2); clf;  q = max(log(P23), l0);
    contourf(flipud(rot90(q)),levs); axis("square"); colormap("jet");
    xlabel("0 < n( i-1 ) < 32 <n( i-1 )>");  ylabel("0 < n(i) < 32 <n(i)>");
    print([fbase,'_lnP23_lev.png'], '-dpng');

figure(3); clf;  q = max(log(P12), l0);
    contourf(flipud(rot90(q)),levs); axis("square"); colormap("jet");
    xlabel("0 < n( i-2 ) < 32 <n( i-2 )>");  ylabel("0 < n( i-1 ) < 32 <n( i-1 )>");
    print([fbase,'_lnP12_lev.png'], '-dpng');

    n1 = ((1:nBins(1)) - 0.5)*dB;
    n2 = ((1:nBins(2)) - 0.5)*dB;
    n3 = ((1:nBins(3)) - 0.5)*dB;
    phi = ((1:nBins(4))-1)*2/nBins(4)-1;

figure(4);  clf; semilogy(n1, P1, '-ob', n2, P2, '-or', n3, P3, '-ok');  set(gca, "fontsize", 16); 
xlabel("n(i-2)/<n(i-2)> (blue),   n(i-1)/<n(i-1)> (red),  n(i)/<n(i)> (black)");  ylabel("probability");
 print([fbase,'_Pn.pdf'], '-dpdf');

figure(5);  clf; plot(phi, P4, '-or'); grid("on"),  set(gca, "fontsize", 16);
xlabel("theta / pi");  ylabel("probability"); axis([-1,1,0.12,0.20])
 print([fbase,'_Ptheta.pdf'], '-dpdf'); 

%return

end %if saveimages

%-- compute mutual information  --

    %-- equilibrium, assuming exponential fit with <n_i> and equidistribution for theta --

    So = 1/log(2) * (1 + log(aa_avg));
    So = [So, log2(2*pi)];

    %-- deviation from equilibrium ---
    
    q =  reshape(P4D, nBins(1)*nBins(2)*nBins(3)*nBins(4), 1);

    ind=find(q>0);
    S4D = - sum(q(ind).*log2(q(ind))) * (dBins(1)*dBins(2)*dBins(3)*dphi)  - sum(So(1:4));
  
    ind=find(P1>0);
    S1 = - sum(P1(ind).*log2(P1(ind))) * dBins(1)  - So(1);

    ind=find(P2>0);
    S2 = - sum(P2(ind).*log2(P2(ind))) * dBins(2) - So(2);

    ind=find(P3>0);
    S3 = - sum(P3(ind).*log2(P3(ind))) * dBins(3) - So(3);

    ind=find(P4>0);
    S4 = - sum(P4(ind).*log2(P4(ind))) * dphi  - So(4);

%   I4D =  S1 + S2 + S3 + log2(2*pi) - S4D;

    I4D =  S1 + S2 + S3 - S4D;

 
    printf('%3d  %7.4f  %7.4f    %7.4f  %7.4f  %7.4f  %7.4f   %s     %4d\n', ...
        i0,  I4D,   S4D,      S1, S2, S3, S4,         fbase,  ntot/1e6);

    clear all

end

%---------------------

