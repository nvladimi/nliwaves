function twomode_P123_slices(fbase_in, fbase_out)
% 
%  reads parameter file for probability distribution (*_post??.param);
%  reads propability distribution P(n1,n2,theta);
%  estimates logarithm of theoretical probability (obsolete);
%  plots averages and slices for numerical and theoretical probability 
%  (comment/uncomment different quantites to show).
%

  global P123;
  global lnP;
  global F;


%--- data ---  

    fname = [fbase_in, '.param'];

    load(fname);  % ('ntot', 'nPhi', 'nN1', 'nN2', 'max_N1', 'max_N2', 'g1', 'g2', 'p1', 'p2' );

    fid = fopen([fbase_in, '.P123'], 'rb');      % P123(nN1, nN2, nPhi+1, 'uint16');

    P123 = fread(fid, nN1*nN2*nPhi, 'uint32');

    fclose(fid);

    P123 = reshape(P123, nN1, nN2, nPhi);

    disp([ntot, sum(sum( sum(P123) )) ]);


%-- setup meshes --

    dn1 = max_N1/nN1; 
    dn2 = max_N2/nN2;
    dphi = 2*pi/nPhi;

    n01 = ((0:nN1-1) + 0.5) * dn1;
    n02 = ((0:nN2-1) + 0.5) * dn2;
    phi0 = ((0:nPhi-1)) * dphi - pi;
   
    [n2, n1, phi] =  meshgrid(n02,n01,phi0);

    %debug: figure(1); q = n1(:,:,1); imagesc(rot90(q)); axis("equal", "off"); colormap("jet");
    %debug: figure(1); plot(1:nPhi, phi(1,1,:), '-or'); grid("on");  set(gca, "fontsize", 16);


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

    disp([alpha, beta]);
    F = F1+F2+F3;
    
%-- to compare with theory, take logarithm of probabiliy and shift origin to zero --

    lnP = log(P123);
    q =  sum( lnP, 3 ) / nPhi;

    q1 = q(1,1);
    q3 = q(2,2);
    q0 = 0.5*(3*q1 - q3);

    lnP = lnP - q0; 

    clear q; clear q0; clear q1; clear q3;

    %debug: plot(n01, q(:,1), '-or', n02, q(1,:), '-ob'); axis([0,0.5,-0.6,0]); grid("on");  set(gca, "fontsize", 16);


%-- data processing --

    % pencils(n01, n02, phi0, fbase_out);
    % prob_theta(phi0, fbase_out)
    images_n1n2(F1(:,:,1), F2(:,:,1), n01, n02, fbase_out)


end

%-----------------------------

function pencils(n01, n02, phi0, fbase_out)

  global lnP;
  global F;

   i1pencils = [0,   0, 10,   0, 10, 20,     0, 10, 20, 30, 40] + 1;
   i2pencils = [0,    5, 0,   10, 5,  0,    20, 15, 10,  5, 0] + 1;      


    N = n01(i1pencils) + 2*n02(i2pencils);
    disp(N);

    m = length(i1pencils);
    Ppencils = zeros(length(phi0), m);
    Fpencils = zeros(length(phi0), m);
    for i=1:m
            q =  lnP( i1pencils(i), i2pencils(i), : );
            qavg = sum(q) / length(q);
            Ppencils(:,i) = q;  %/qavg - 1;
            q =  F( i1pencils(i), i2pencils(i), : );
            qavg = sum(q) / length(q);
            Fpencils(:,i) = q; %/qavg - 1;
    end;


    ip0 = phi0/pi;

    figure(1); clf
    
      plot(...
          ip0, Ppencils(:,1),      '--ok', ...
          ip0, Ppencils(:,2:3),    '--or', ...
          ip0, Ppencils(:,4:6),    '--ob', ...
          ip0, Ppencils(:,7:11),   '--om', ...
          ip0, Fpencils(:,1),      '-k', ...
          ip0, Fpencils(:,2:3),    '-r', ...
          ip0, Fpencils(:,4:6),    '-b', ...
          ip0, Fpencils(:,7:11),   '-m');

     %axis([-1,1, -0.1,0.1]);
     set(gca, "fontsize", 16);  title(fbase_out);
     print([fbase_out,'_phi_pensils.pdf'], '-dpdf');

end


%-----------------------------

function prob_theta(phi0, fbase_out)

   global P123;

   P3 = squeeze(sum(sum(P123,2),1));

   P3avg = sum(P3)/length(P3);

   P3  = P3/P3avg - 1;
   ip0 = phi0/pi;

   figure(1); clf
   plot(ip0, P3, '-ok');
   grid("on"); title(fbase_out);
   %axis([-1,1, -0.1,0.1]);
   set(gca, "fontsize", 16);
   print([fbase_out,'_probtheta.pdf'], '-dpdf');

end
%-------------------------------

function images_n1n2(f1, f2, n01, n02, fbase_out)

    global lnP;

    s = size(lnP);

    lnP12 =  sum( lnP, 3 ) / s(3);

    n1max = 80; n2max =40;

    lnP12 = lnP12(1:n1max,1:n2max);
    f1    = f1(1:n1max,1:n2max);
    f2    = f2(1:n1max,1:n2max);
    n01   = n01(1:n1max);
    n02   = n02(1:n2max);

    levs =(-10:1:0);

    % figure(1); clf;
    % q = max(lnP12, -10);
    % contourf(flipud(rot90(q)),levs); axis("equal", "off"); colormap("jet");
    % print([fbase_out,'_lnP_lev.png'], '-dpng');

    % figure(2); clf;
    % q = max(f1, -10);
    % contourf(flipud(rot90(q)),levs); axis("equal", "off"); colormap("jet");
    % print([fbase_out,'_f1_lev.png'], '-dpng');

    levs =(-0.3:0.01:0.3);
    
    % figure(1);
    % q = max(lnP12-f1, -0.5);
    % contourf(flipud(rot90(q)),levs); axis("equal", "off"); colormap("jet");
    % print([fbase_out,'_lev.png'], '-dpng');
    %
    %  figure(2);
    %  q = max(f2, -10);
    %  contourf(flipud(rot90(q)),levs); axis("equal", "off"); colormap("jet");caxis([-0.2,0.2]);
    %  print([fbase_out,'_f2_lev.png'], '-dpng');


    % figure(1); clf
    % q = lnP12 - f1;
    % m = max(max(q)), cmap=[-m,m]; 
    % imagesc(rot90(q), cmap); axis("equal", "off"); colormap("jet");
    % print([fbase_out,'_lnPdiff_img.png'], '-dpng');

    % figure(2); clf;
    % q = f2;
    % m = max( abs(min(min(q))), max(max(q)) ), cmap =[-m,m];
    % imagesc(rot90(q), cmap); axis("equal", "off"); colormap("jet");
    % print([fbase_out,'_f2_img.png'], '-dpng');


    figure(1); clf;
    plot(n01, lnP12(:,1), '--or' , n01, f1(:,1), '-r',  ...
	 n02, lnP12(1,:), '--ob',  n02, f1(1,:), '-b',  ...
	 n02, lnP12(20,:), '--om',  n02, f1(20,:), '-m' );
    grid("on"); title(fbase_out);  set(gca, "fontsize", 16); 
    % print([fbase_out,'_lnP.pdf'], '-dpdf');


    figure(2); clf
    plot(n01, lnP12(:,1)-f1(:,1), '--or' , n01, f2(:,1), '-r', ...
	 n02, lnP12(1,:)-f1(1,:), '--ob', n02, f2(1,:), '-b',  ...
	 n02, lnP12(20,:)-f1(20,:), '--om', n02, f2(20,:), '-m' );
grid("on"); title(fbase_out);  set(gca, "fontsize", 16); axis([0,8,-0.4,0.4])
     print([fbase_out,'_lnPdiff.pdf'], '-dpdf');




end
%-------------------------------

