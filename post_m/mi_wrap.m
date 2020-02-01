function mi_wrap(fbase, fbaseout, fnums, N, N0, ndiag, amp0, nAmp, nPhi, m1, m2)
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

%--------------------------------

n = length(fnums)*ndiag;  % total number of sumples

if (m1 >= 0) im1 = m1+1;  else im1 = N + m1 + 1; end      % mode index
if (m2 >= 0) im2 = m2+1;  else im2 = N + m2 + 1; end      % mode index

[P12amp, Pdphi] = mi_core(fbase, fnums, N, N0, ndiag, amp0, nAmp, nPhi, im1, im2);

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

end
