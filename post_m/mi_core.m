function mi_core(fbase, fbaseout, fnums, N, N0, ndiag, amp0, nAmp, nPhi,  m1, m2);
%
% fbase       string, base for input files;
% fbaseout    string, base for output files;
% fnums       array of files numbers to process;
% N           saved number of modes in input file; 
% N0          total number of modes in original simulation;
% ndiag       number of samples per file;
% amp0        max amplitude of a mode;
% nAmp        number of bins for amplitudes;
% nPhi        number of bins for phase (must be even!)
% m1          mode, kx;
% m2          mode, ky;
%


iunit=i;

n = length(fnums)*ndiag;  % total number of sumples

if (m1 >= 0) im1 = m1+1;  else im1 = N + m1 + 1; end      % mode index
if (m2 >= 0) im2 = m2+1;  else im2 = N + m2 + 1; end      % mode index


%###########################################################






%###########################################################


%-- create empty array to store probabilities --

disp('-- create probability array --');
tic

  nA = nAmp;
  nF = nPhi + 1;

  nbins_A = nAmp - 1;
  nbins_phi = round( nPhi/2 );
  NN  = N*N;

  P12 = zeros(nA, nF, nA*nF*NN, 'uint16');


  % dk=2*pi/L;
  k = 0:N-1;
  k0 = floor((N-1)/2); 
  k = circshift(k-k0,[0,-k0]);
  [Kx, Ky] = meshgrid(k,k);
  M = (-1).^(Kx + Ky);
  kk = Kx.*Kx + Ky.*Ky;

  Ampk = ones(N,N) * amp0 ./ kk;    % amplitutude for conversion into index

toc
  
%-------------------------------------------------------------------------------------
%-- read files, compute amplitude and phase, --
%-- convert amplitude and phase into bin indexes and store in propability arrays  --

  nkavg = zeros(N,N);
  dfavg = zeros(N,N);
  aaavg = zeros(N,N);

  
for fnum=fnums

  disp('-- read in data -- ');

  tic

  fname = [fbase, '.klo.',  num2str(fnum,'%04d')]

  fid=fopen(fname, 'rb');
  f=fread(fid, N*N*2*ndiag, 'double'); 
  fclose(fid);

  f = reshape(f,2,N,N,ndiag) * N0*0.5; 
  f = f(1,:,:,:) - iunit*f(2,:,:,:); 
  f = reshape(f,N,N,ndiag);  

  for i=1:ndiag
      f(:,:,i) = f(:,:,i) .* M;
  end;

  toc


  disp('-- add to probapility arrays --');
  tic

  iNN = transpose(1:NN);

  aaavg1=zeros(N,N);
  dfavg1=zeros(N,N);
  nkavg1=zeros(N,N);

  for i=1:ndiag

          f1 = f(:,:,i);
          p1 = angle(f1);
          a1 = abs(f1);

          f0 = conj(f1(im1,im2));
          p0 = p1(im1,im2);
          dp = p1 - p0;

          ind = find (dp >  pi);  dp(ind) = dp(ind) - 2*pi;
          ind = find (dp < -pi);  dp(ind) = dp(ind) + 2*pi;

          nkavg1 = nkavg1 + a1.*a1;
          aaavg1 = aaavg1 + f1*f0;
          dfavg1 = dfavg1 + dp;


	  iA = round( a1./Ampk * nbins_A) + 1;
	  iF = round( p1./pi * nbins_phi) + nbins_phi + 1;

          iA  = min(iA, nA);

	  mA   = iA(im1,im2);
	  mF   = iF(im1,im2);

          iA  = reshape(iA, [NN,1]);
          iF  = reshape(iF, [NN,1]);

          ind=sub2ind([nA, nF, NN], iA, iF, iNN);

          P12(mA, mF, ind) = P12(mA, mF, ind) + uint16(1);
  end

  nkavg = nkavg + nkavg1/n;
  aaavg = aaavg + aaavg1/n;
  dfavg = dfavg + dfavg1/n;

  toc

end % cycle over files
clear f

%------------------------------------------------------------

disp('-- probabilities for phase difference and amplitudes --');
tic

   P12 = reshape(P12, nA, nF, nA, nF, N, N);

   P12amp  = squeeze(sum(sum(P12,2),4));  %<--  P12amp(nA, nA, N, N);

   P12phi  = squeeze(sum(sum(P12,1),3));  %<--  P12phi(nF, nF, N, N);


   ind=zeros(nF-1,nF-1);
   for i=1:nF-1
       ind = ind + i*circshift(eye(nF-1), [i-1,0]);
   end
   ind(:,end+1)=ind(:,1);
   ind(end+1,:)=ind(1,:);

   Pdphi  = zeros(nF-1, N, N);
 
   for i1=1:nF
      for i2=1:nF
	    i = ind(i1,i2);
            q = squeeze(P12phi(i1,i2,:,:)) + squeeze(Pdphi(i,:,:));
            Pdphi(i,:,:)  =  q;
      end
   end

toc

%##############################################################

      
[nA,nA,N,N] = size(P12amp); 
[nF,   N,N] = size(Pdphi); 
  
NN = N*N;

fnameout = [fbaseout, '.param'];

save( fnameout, 'fbase', 'm1', 'm2', 'n', 'N', 'nA', 'nF', 'amp0', 'nAmp', 'nPhi', 'ndiag', 'fnums' ); 

fid = fopen([fbaseout, '.P12amp'], 'wb');
fwrite(fid, P12amp, 'uint16');
fclose(fid);

fid = fopen([fbaseout, '.Pdphi'], 'wb');
fwrite(fid, Pdphi, 'uint32');
fclose(fid);

fid = fopen([fbaseout, '.nkavg'], 'wb');
fwrite(fid, nkavg, 'double');
fclose(fid);

fid = fopen([fbaseout, '.aamavg'], 'wb');
fwrite(fid, abs(aaavg), 'double');
fclose(fid);

fid = fopen([fbaseout, '.aapavg'], 'wb');
fwrite(fid, angle(aaavg), 'double');
fclose(fid);


fid = fopen([fbaseout, '.dfavg'], 'wb');
fwrite(fid, dfavg, 'double');
fclose(fid);


%-------------------------------

disp('-- compute mutual information for amplitudes --');
tic;

   P1amp  = squeeze(sum(P12amp,1));

   P12amp = reshape(P12amp,[nA*nA, NN]);
   P1amp  = reshape(P1amp, [nA, NN]);

   S12 = zeros(NN,1);
   S1  = zeros(NN,1);

   for i=1:NN
	   q=squeeze(P12amp(:,i));
           ind=find(q>0);
           S12(i) = sum(q(ind).*log(q(ind)));
   end;

   for i=1:NN
	   q=squeeze(P1amp(:,i));
           ind=find(q>0);
           S1(i) = sum(q(ind).*log(q(ind)));
   end;

   S12 = S12/n - log(n);
   S1  = S1/n  - log(n);

   S12 = reshape(S12,[N,N]);
   S1  = reshape(S1, [N,N]);

   S2 = S1(im1,im2);

   I12 = S12 - S1 - S2; 

   P12amp = reshape(P12amp,[nA,nA,N,N]);
 
toc

     
disp('-- compute average amplitude and phase and their deviation --');
tic

 
  Pamp = reshape(sum(P12amp,1), [nA, NN]);
  Pphi = reshape(Pdphi, [nF, NN]);

  iiA = transpose(1:nA);
  iiF = transpose(1:nF);

  iAavg=zeros(NN,1);
  iAdev=zeros(NN,1);

  iFavg=zeros(NN,1);
  iFdev=zeros(NN,1);


  for i=1:NN

	  c    = Pamp(:,i);
          csum = squeeze(sum(c));
          iavg = sum(c.*iiA)/csum;
          idev = sum((c.*(iiA - iavg).^2))/csum; 
          iAavg(i) = iavg;
          iAdev(i) = sqrt(idev);

	  c    = Pphi(:,i);
          csum = sum(c);
          iavg = sum(c.*iiF)/csum;
          idev = sum(c.*(iiF-iavg).^2)/csum; 
          iFavg(i) = iavg;
          iFdev(i) = sqrt(idev);

  end 

  iAavg = reshape(iAavg, [N,N]); 
  iAdev = reshape(iAdev, [N,N]); 

  iFavg = reshape(iFavg, [N,N]); 
  iFdev = reshape(iFdev, [N,N]); 

toc


fnameout = [fbaseout, '.mi'];

save( fnameout, 'N', 'I12', 'S12', 'S1', 'S2', 'iAavg', 'iAdev', 'iFavg', 'iFdev' ); 

      
return


end
