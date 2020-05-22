function twomode_P123_slices(fbase_in, fbase_out)


n1pencils = [0,   0, 10,   0, 10, 20,     0, 10, 20, 30, 40] + 1;
n2pencils = [0,    5, 0,   10, 5,  0,    20, 15, 10,  5, 0] + 1;      

phi_slices = [1, 8, 16, 24];

%--- data ---  

   fname = [fbase_in, '.param'];

    load(fname);  % ('ntot', 'nPhi', 'nN1', 'nN2', 'max_N1', 'max_N2', 'g1', 'g2', 'p1', 'p2' );

    fid = fopen([fbase_in, '.P123'], 'rb');      % P123(nN1, nN2, nPhi+1, 'uint16');

    P123 = fread(fid, nN1*nN2*nPhi, 'uint16');

    fclose(fid);

    P123 = reshape(P123, nN1, nN2, nPhi);

    disp([ntot, sum(sum( sum(P123) )) ]);
	 
%-- theory --

    dn1 = max_N1/nN1; 
    dn2 = max_N2/nN2;
    dphi = 2*pi/nPhi;

    n01 = ((0:nN1-1) + 0.5) * dn1;
    n02 = ((0:nN2-1) + 0.5) * dn2;
    phi0 = ((0:nPhi-1)) * dphi - pi;
   
    [n2, n1] =  meshgrid(n01,n02);

    P12 = sum( P123, 3 );
    P12 = P12 / sum(sum(P12)) / (dn1*dn2);

    N12 = n1 + 2*n2;
    %N12 = N12 / sum(sum(N12)) / (dn1*dn2);

%-- images --

p=log(P12);
p = -N12;

p = p(1:100,1:50);
p = max(p, -10);

minP = min(min(p))
maxP = max(max(p))
cmap = [-10,0]; levs =(-10:1:0);

figure(1);
imagesc(rot90(p), cmap); axis("equal", "off"); colormap("jet");
print([fbase_out,'_img.png'], '-dpng');


figure(2);
contourf(rot90(p),levs); axis("equal", "off"); colormap("jet");
print([fbase_out,'_lev.png'], '-dpng');




return






%-- pencils --

    P3  = squeeze(sum( sum (P123, 2), 1 ));
    P3 = P3 / sum(P3) / dphi;
    P3avg = sum(P3)/length(P3);
    dP3 = P3/P3avg - 1;


    m = length(n1pencils);
    Ppencils = zeros(nPhi, m);
    for i=1:m
      Ppencils(:,i) = P123(n1pencils(i), n2pencils(i),:);
    end;

    Ppencils_avg = sum(Ppencils,1)/nPhi;
    Ppencils_avg  = repmat(Ppencils_avg, [nPhi, 1]);
    dPpencils = Ppencils./Ppencils_avg - 1;

    N = n01(n1pencils) + 2*n02(n2pencils);


    figure(1);
    ip0=phi0/pi;
    plot(ip0, dP3, '-ok'); grid("on");
    set(gca, "fontsize", 16);
    %print([fbase_out,'_phi_avg.pdf'], '-dpdf');


    figure(2);
    ip0=phi0/pi;
    plot(ip0, dP3, '-ok',  ...
         ip0, dPpencils(:,1),      '-r', ...
         ip0, dPpencils(:,2:3),    '-b', ...
         ip0, dPpencils(:,4:6),    '-m', ...
         ip0, dPpencils(:,7:11),   '-c'); axis([-1,1, -0.5,0.5]);
     set(gca, "fontsize", 16);
     %print([fbase_out,'_phi_pensils.pdf'], '-dpdf');

%-- end pencils ---

return






figure(3);
x=phi0/pi; y1=P12(5,:); y2=P12(10,:); y3=P12(20,:); y4=P12(25,:);
plot(x, y1, '-ok',  x,y2, '-ob',  x,y3, '-or',  x,y4, '-o', 'color', [0 0.5 0]);
set(gca, "fontsize", 20);
%print([fbase_out,'_plot.pdf'], '-dpdf');

return
clf();

figure(1);  cmap=[-1,1]*1.15*maxK;
imagesc(rot90(K),cmap); axis("equal", "off"); colormap("jet"); %set(gca, "fontsize", 20);
print([fbase_out,'_img.png'], '-dpng');

figure(2);  levs = (-1:0.1:1)*maxK;
contourf(rot90(K),levs); axis("equal", "off"); colormap("jet"); set(gca, "fontsize", 20);


figure(3);
x=phi0/pi; y1=K(5,:); y2=K(10,:); y3=K(20,:); y4=K(25,:);
plot(x, y1, '-ok',  x,y2, '-ob',  x,y3, '-or',  x,y4, '-o', 'color', [0 0.5 0]);
set(gca, "fontsize", 20);
print([fbase_out,'_plot.pdf'], '-dpdf');

end

%---------------------

