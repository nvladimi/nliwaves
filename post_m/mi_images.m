function mi_images

%fbaseout = '../mi/a5_m06m12_post01';
%fbaseout = '../mi/a5_m06m12_post02';
%fbaseout = '../mi/a5_m06m12_post03';


%fbaseout = '../mi/a5_m06m12_post03';

%fbaseout = '../mi/a5_m11m13_post03';

fbaseout = '../mi/a1_m06m10_post02';



load([fbaseout, '.param']);
load([fbaseout, '.mi']);


%------------------------------


im1 = m1+1;               % mode index
im2 = m2+1;               % mode index


fid = fopen([fbaseout, '.P12amp'], 'rb');
P12amp = fread(fid, nA*nA*N*N, 'uint16');
fclose(fid);

fid = fopen([fbaseout, '.Pdphi'], 'rb');
Pdphi = fread(fid, nF*N*N, 'uint16');
fclose(fid);

%-------------------------------

   P12amp = reshape( P12amp, [nA,nA, N,N] ); 
   Pdphi  = reshape( Pdphi,  [nF,    N,N] ); 

   P1amp  = squeeze(sum(P12amp,1));

q=iFavg; disp([ min(min(q)), max(max(q)), sum(sum(q))/N/N ])
figure(1); imagesc(rot90(fftshift(iFavg)), [8.45,8.55]); axis("square","off"); colormap("jet"); 

q=iAavg; disp([ min(min(q)), max(max(q)), sum(sum(q))/N/N ])
figure(2); imagesc(rot90(fftshift(iAavg)), [10,20]); axis("square","off"); colormap("jet"); 

% q=S12; disp([ min(min(q)), max(max(q)), sum(sum(q))/N/N ])
% figure(3); imagesc(rot90(fftshift(S12)), [-7,-5]); axis("square","off"); colormap("jet"); 

% q=S1; disp([ min(min(q)), max(max(q)), sum(sum(q))/N/N ])
% figure(4); imagesc(rot90(fftshift(S1)), [-4,-3]); axis("square","off"); colormap("jet"); 

q=I12; disp([ min(min(q)), max(max(q)), sum(sum(q))/N/N ])
figure(5); imagesc(rot90(fftshift(I12)), [0,0.02]); axis("square","off"); colormap("jet"); 



%p = squeeze( P12amp(:,:, im1,im2) ); 
%figure(2); imagesc(rot90(fftshift(p))); axis("square","off"); colormap("jet"); 

%p = squeeze( P12amp(:,:, 5,5) ); 
%figure(3); imagesc(rot90(log(p)); axis("square","off"); colormap("jet"); 





%  clf(); imagesc(fftshift(iAavg)); axis("square", "off"); colormap("jet");  print('IMG/iAavg.png', '-dpng');
%  clf(); imagesc(fftshift(iAdev)); axis("square", "off"); colormap("jet");  print('IMG/iAdev.png', '-dpng');

%  clf(); imagesc(fftshift(iFavg)); axis("square", "off"); colormap("jet");  print('IMG/iFavg.png', '-dpng');
%  clf(); imagesc(fftshift(iFdev)); axis("square", "off"); colormap("jet");  print('IMG/iFdev.png', '-dpng');

%  clf(); imagesc(fftshift(I12)); axis("square", "off"); colormap("jet");  print('IMG/I12.png', '-dpng');
%  clf(); imagesc(fftshift(S12)); axis("square", "off"); colormap("jet");  print('IMG/S12.png', '-dpng');
%  clf(); imagesc(fftshift(S1));  axis("square", "off"); colormap("jet");  print('IMG/S1.png', '-dpng');


%  imagesc(fftshift(I12), [0,1]); axis("square","off"); colormap("jet"); 

%plot(1:16, Pdphi(:,2,4), '-o'); set(gca, "fontsize", 20);


end
