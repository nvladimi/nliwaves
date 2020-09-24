function fibo_mi3D(fbase, saveimages)
%
%  reads parameter file for probability distribution
%  reads propability distribution P3D(n2,n3,theta)
%  computes MI and entropies for this distributon
% 
%  header
%  #1.fname  2.samples(M)  3.S3D    4.S2     5.S3     6.S4    7.I3D
%


%--- data ---  

    fname = [fbase, '.param'];
    load(fname);  % ('i0', 'nBins', 'dB', 'ntot',  'aa_avg'); 

    fname = [fbase, '.P3D'];
    fid = fopen(fname, 'rb');      % P3D(nN2, nN3, nPhi, 'uint32');

    P3D = fread(fid, nBins(2)*nBins(3)*nBins(4), 'uint32');

    fclose(fid);
 
    P3D = reshape(P3D, nBins(2), nBins(3), nBins(4));


    P2  = squeeze(sum( sum(P3D, 3), 2));
    P3  = squeeze(sum( sum(P3D, 3), 1));
    P4  = squeeze(sum( sum(P3D, 2), 1));

    P23 = squeeze( sum(P3D,3) );

	 
%-- normalise distributions to one  --

    dBins = aa_avg*dB;

    dphi = 2*pi/nBins(4);

    P3D  = P3D/ntot/(dBins(2)*dBins(3)*dphi);
    P2   = P2/ntot/dBins(2);
    P3   = P3/ntot/dBins(3);
    P4   = P4/ntot/dphi;

    P23  = P23/ntot/dBins(2)/dBins(3);

if saveimages

    warning('off');

    %l0=-5; levs =(l0:1:l0+15);
    %l0=-15; levs =(l0:1:10+15);
    l0=-12; levs =(l0:1:5);


    % figure(1);  imagesc(log(rot90(P23))); axis("square", "off"); colormap("jet");
 
    figure(1); clf;  q = max(log(P23), l0);
    contourf(flipud(rot90(q)),levs); axis("square"); colormap("jet");
    xlabel("0 < n( i-1 ) < 32 <n( i-1 )>");  ylabel("0 < n(i) < 32 <n(i)>");
    %print([fbase,'_lnP23_lev.png'], '-dpng');

    n2 = ((1:nBins(2)) - 0.5)*dB;
    n3 = ((1:nBins(3)) - 0.5)*dB;
    phi = ((1:nBins(4))-1)*2/nBins(4)-1;

    figure(2);  clf; semilogy(n2, P2, '-or', n3, P3, '-ok');  set(gca, "fontsize", 16); 
    xlabel("n(i-1)/<n(i-1)> (red),  n(i)/<n(i)> (black)");  ylabel("probability");
    %print([fbase,'_Pn.pdf'], '-dpdf');

    figure(3);  clf; plot(phi, P4, '-or'); grid("on"),  set(gca, "fontsize", 16);
    xlabel("phi / pi");  ylabel("probability"); axis([-1,1,0.12,0.20])
    print([fbase,'_Pphi.pdf'], '-dpdf'); 

end %if saveimages

%-- compute mutual information  --

    %-- equilibrium, assuming exponential fit with <n_i> and equidistribution for theta --

    So = 1/log(2) * (1 + log(aa_avg));
    So = [So, log2(2*pi)];

    %-- deviation from equilibrium ---
    
    q =  reshape(P3D, nBins(2)*nBins(3)*nBins(4), 1);

    ind=find(q>0);
    S3D = - sum(q(ind).*log2(q(ind))) * (dBins(2)*dBins(3)*dphi)  - sum(So(2:4));
  
    ind=find(P2>0);
    S2 = - sum(P2(ind).*log2(P2(ind))) * dBins(2) - So(2);

    ind=find(P3>0);
    S3 = - sum(P3(ind).*log2(P3(ind))) * dBins(3) - So(3);

    ind=find(P4>0);
    S4 = - sum(P4(ind).*log2(P4(ind))) * dphi  - So(4);

%   I3D =  S2 + S3 + log2(2*pi) - S3D;

    I3D =  S2 + S3 - S3D;

 
    printf('%s  %4d    %7.4f     %7.4f  %7.4f %7.4f   %7.4f\n', ...
	   fbase,  ntot/1e6,   S3D,   S2, S3, S4,  I3D);
 
end

%---------------------

