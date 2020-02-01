function [P12amp, Pdphi] = mi_core(fbase, fnums, N, N0, n, amp0, nAmp, nPhi,  m1, m2);
  

iunit=i;


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

%-- read files, compute amplitude and phase, --
%-- convert amplitude and phase into bin indexes and store in propability arrays  --

  
for fnum=fnums

  disp('-- read in data -- ');

  tic

  fname = [fbase, '.klo.',  num2str(fnum,'%04d')]

  fid=fopen(fname, 'rb');
  f=fread(fid, N*N*2*n, 'double'); 
  fclose(fid);


  f = reshape(f,2,N,N,n); 
  f = f(1,:,:,:) - iunit*f(2,:,:,:); 
  f = reshape(f,N,N,n);  

  for i=1:n
      f(:,:,i) = f(:,:,i) .* M;
  end;

  amp = abs(f)*N0*0.5;
  phi = angle(f);
  
  toc


  disp('-- add to probapility arrays --');
  tic

  iNN = transpose(1:NN);

  for i=1:n

	  iA = round( amp(:,:,i)./Ampk * nbins_A) + 1;
	  iF = round( phi(:,:,i)./pi * nbins_phi) + nbins_phi + 1;

          iA  = min(iA, nA);

	  mA   = iA(m1,m2);
	  mF   = iF(m1,m2);

          iA  = reshape(iA, [NN,1]);
          iF  = reshape(iF, [NN,1]);

          ind=sub2ind([nA, nF, NN], iA, iF, iNN);

          P12(mA, mF, ind) = P12(mA, mF, ind) + uint16(1);
  end

  toc

  P12 = reshape(P12, nA, nF, nA, nF, N, N);

end % cycle over files

clear n
clear f

disp('-- probabilities for phase difference and amplitudes --');
tic

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
      
return

%-----------  UNUSED ----------------------

disp('-- compute mutial information in 4D --');
tic;

   P12 = reshape(P12,[nA*nF*nA*nF, NN]);
   P1  = reshape(sum(sum(P12,1),2), [nA*nF, NN]);  %<-- verify!!!

   s12 = zeros(NN,1);
   s1  = zeros(NN,1);

   for i=1:NN
	   q=squeeze(P12(:,i));
           ind=find(q>0);
           s12(i) = sum(q(ind).*log(q(ind)));
   end;

   for i=1:NN
	   q=squeeze(P1(:,i));
           ind=find(q>0);
           s1(i) = sum(q(ind).*log(q(ind)));
   end;

   s12 = s12/ntot - log(ntot);
   s1  = s1/ntot  - log(ntot);

   s12 = reshape(s12,[N,N]);
   s1  = reshape(s1, [N,N]);

   s2 = s1(m1,m2);

   I12 = s12 - s1 - s2; 

toc


end
