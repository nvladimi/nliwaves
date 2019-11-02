clear all

n=2000; 
N=8;                     % saved number of modes 
N0=256;                  % total number of modes
fname='n0400/n0400.klo'; 
tmax=2;

iunit=i;

%-- read file of psi_k for low k --

fid=fopen(fname, 'rb');
f=fread(fid,N*N*2*n, 'double'); 
fclose(fid); 

f=reshape(f,2,N,N,n); 
f=f(1,:,:,:) + iunit*f(2,:,:,:); 
f=reshape(f,N,N,n);  

h=abs(f)*N0;                 %-- amplitude of modes vs time
fh=fft(h,n,3)/n;             %-- fourier in time
fh=abs(fh);
 

%-- find condensate level and Bogolubov frequency --

A0 = fh(1,1,1)
Omega = A0*sqrt(2);
T=2*pi/Omega;
t=(1:n)/n*tmax;  
dw=2*pi/tmax/Omega;
w=(0:n-1)*dw;

%-- scale time-spectra with condensate --

fh = fh/A0;


%-- wipe out lowest five points, non-periodicity of time domain -- 

fh(:,:,1:6)=0; fh(:,:,end-4:end)=0;

%-- store condensate time-spectrum --

h0 = reshape(fh(1,1,:),1,n); 

%-- find peaks in spectra of individual modes --

for i=1:N
  for j=1:N
     q=reshape(fh(i,j,:), 1,n);
     [qm,im] = max(q);
     fhmax(i,j) = qm;
     wmax(i,j) = im;
  end
end


wmax=fftshift(wmax); wmax=wmax(2:end,2:end);
fhmax=fftshift(fhmax); fhmax=fhmax(2:end,2:end);

wmax = wmax*dw;

Frequencies = wmax(3:5,3:5)
Amplitudes = fhmax(3:5,3:5)


%figure(1, 'Name',  'Frequency'); imagesc(wmax); axis('equal', 'tics', 'off')
%figure(2, 'Name',  'Amplitude'); imagesc(fhmax, [0,0.35]); axis('equal', 'tics', 'off')

%-- plot first two modes vs time --

h1=reshape(h(1,2,:),1,n); 
h2=reshape(h(2,2,:),1,n);

figure(3, 'Name', 'modes k=(0,1) and k=(1,1)'); 
plot(t,h1/A0, t,h2/A0); %axis([0,0.5])
%plot(t/T, h2/A0); axis([0,2])

%-- plot spectra of individual modes --

h1=reshape(fh(1,2,:),1,n); 
h2=reshape(fh(2,1,:),1,n);  
h3=reshape(fh(end,1,:),1,n);   
h4=reshape(fh(1,end,:),1,n); 

#figure(4, 'Name', 'spectrum of modes |k| = 1' ); 
#plot(w,h1,w,h2,w,h3,w,h4); axis([0,10])


h1=reshape(fh(2,2,:),1,n); 
h2=reshape(fh(2,end,:),1,n);  
h3=reshape(fh(end,2,:),1,n);   
h4=reshape(fh(end,end,:),1,n); 

%figure(5, 'Name', 'spectrum of modes |k| = sqrt(2)' ); 
%plot(w,h1,w,h2,w,h3,w,h4); axis([0,10])

%figure(6, 'Name', 'spectrum of condensate');
%plot(w, h0, 'o-'); axis([0,10])

%---------------------

[hm1,im1]=max(h0);
printf('Consentate: w=%f,   A/A0=%e\n', im1*dw, hm1);

h0(im1)=0; h0(end-im1+2)=0;

[hm2,im2]=max(h0);

printf('Consentate: w=%f,   A/A0=%e\n', im2*dw, hm2);


%--------------------

%-- print data for "summary" file --

%_1.N  2.A0  3.(-1,-1)  4.(-1,0)  5.(-1,1)  6.(0,-1)  7.(0,0)  8.(0,1)  9.(1,-1)  10.(1,0)  11.(1,1)

disp([ A0, reshape(Amplitudes,1,9) ]);

%_1.N  2.A0    3.w1    4.w2      5.wmax1   6.Amax1/A0    7.wmax2   8.Amax2/A0


printf("%6.3f   %6.4f  %6.4f   %8.6f  %12.6e   %8.6f  %12.6e\n", ...
         A0,  Frequencies(1,2),  Frequencies(1,3), im1*dw, hm1,  im2*dw, hm2)
