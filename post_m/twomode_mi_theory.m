function twomode_mi_theory(fbase)
%
%  reads parameter file for probability distribution (*_post??.param)
%  creates propability distribution based on theoretical assumptions
%  computes MI and entropies for generated distributon
%
 
%--- data ---  

    fname = ['POST/post02/', fbase, '_post02.param'];

    load(fname);  % ('ntot', 'nPhi', 'nN1', 'nN2', 'max_N1', 'max_N2', 'g1', 'g2', 'p1', 'p2' );

%-- playground: what if parameters were different? ---

% max_N1 = max_N1*2;
% max_N2 = max_N2*2;

% nN1 = nN1*2;
% nN2 = nN2*2;

% g1 = 0.06;
% g2 = 0.06;

%-- setup meshes --

    dn1 = max_N1/nN1; 
    dn2 = max_N2/nN2;
    dphi = 2*pi/nPhi;

    n01 = ((0:nN1-1) + 0.5) * dn1;
    n02 = ((0:nN2-1) + 0.5) * dn2;
    phi0 = ((0:nPhi-1)) * dphi - pi;
   
    [n2, n1, phi] =  meshgrid(n02,n01,phi0);

%-- theory ---

    T1 = 0.5 * p1 / g1;
    T2 = p2 / g2;
    dT = 2 * (T1 - T2);
    T  = 0.5 * (p1 + 2*p2) / (g1 + g2);
    chi = (g1 + g2)^3 / (p1 + 2*p2);

    c1  = 1 + 2*p2/p1;
    c2 = g1*g2 / (g1 + g2)^2;
    c3 = (g1 + g2) / (2*g1 + g2);
    c4 = (g1 + g2) / (4*g1);


    if (chi < 1)
         alpha  = 8*c1*c2 * dT/T;
         beta   = c2/8;
         F1 = -(n1 + 2*n2);
         F2 =  alpha * ( n1.*n2 - beta*n1.*n1);
         F3 = -alpha * sqrt(chi/2) * n1.*sqrt(n2) .* sin(phi); 
    else
         alpha  = c3/chi * T*dT/(T1*T2);
         beta = c4;
         F1 = -(n1*T/T1 + 2*n2*T/T2);
         F2 =  alpha * ( n1.*n2 - beta*n1.*n1 );
         F3 =  alpha * sqrt(2*chi) * n1.*sqrt(n2) .* sin(phi); 
    end;

    F = F1 + F2 + F3;
    P123 = exp(F);
    n123 = sum(sum(sum(P123)));

    P1  = squeeze(sum( sum(P123,3), 2 ));
    P2  = squeeze(sum( sum(P123,3), 1 ));
    P3  = squeeze(sum( sum(P123,2), 1 ));

%-- normalise distributions to one  --

    dn1 = max_N1/nN1; 
    dn2 = max_N2/nN2;
    dphi = 2*pi/nPhi;

    P123 = P123/n123/(dn1*dn2*dphi);
    P1   = P1/n123/dn1;
    P2   = P2/n123/dn2;
    P3   = P3/n123/dphi;

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

    I12 =  S1 + S2  + log2(2*phi) - S123;


%-- print parameters and theoretical entropies --

   printf("%s | %4.2f  %4.2f  %5.2e  %5.2e | ", fbase,  g1, g2, p1, p2);
   printf("%8.2e  %8.3f  %6.3f  %8.5f  %8.5f  | ", ntot, chi, dT/T,  alpha,  beta);
   printf('%7.4f  %7.4f  %7.4f  %7.4f  %9.4e\n', S123, S1, S2, S3,  I12);
 
end

%-----------------------------
