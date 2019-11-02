function guide_pearson

N=8;                     % saved number of modes 
N0=256;                  % total number of modes

%fname='n0400a.klo';  A0=19.75; n=4000; tmax=0.4;
 fname='n3600a.klo';  A0=59.73; n=4000; tmax=0.2;

iunit=i;

%T=2*pi/sqrt(2)/A0
%Omega = 2*pi/T
tdata=(1:n)/n*tmax;  


kx=0; ky=1;  i1=1; i2=1; j1=2; j2=8;  fnameext = 'k01';
%kx=1; ky=0;  i1=2; i2=8; j1=1; j2=1;  fnameext = 'k10';
%kx=1; ky=1;  i1=2; i2=8; j1=2; j2=8;  fnameext = 'k11';
%kx=1; ky=-1; i1=2; i2=8; j1=8; j2=2;   fnameext = 'k-11';

%kx=3; ky=3;  i1=4; i2=6; j1=4; j2=6;  fnameext = 'k33';
%kx=3; ky=-3; i1=4; i2=6; j1=6; j2=4;  fnameext = 'k-33';
%kx=0; ky=3;  i1=1; i2=1; j1=4; j2=6;  fnameext = 'k03'; 
%kx=3; ky=0;  i1=4; i2=6; j1=1; j2=1;  fnameext = 'k30';
		  
kk = kx^2 + ky^2
  
Omega = sqrt(2*A0^2*kk + kk^2); T=2*pi/Omega


%-- read file of psi_k for low k --

fid=fopen(fname, 'rb');
f=fread(fid,N*N*2*n, 'double'); 
fclose(fid); 

f=reshape(f,2,N,N,n); 
f=f(1,:,:,:) + iunit*f(2,:,:,:); 
f=reshape(f,N,N,n);  

h=abs(f)*N0;     %-- amplitude of modes vs time
a=arg(f);

theta = reshape( 2*a(1,1,:) - a(i1,j1,:) - a(i2,j2,:), n);
amplitude = 0.5*reshape( h(i1,j1,:) + h(i2,j2,:), n);

ind = find(theta>2*pi);
theta(ind)=theta(ind)-2*pi;
ind = find(theta<-2*pi);
theta(ind)=theta(ind)+2*pi;
ind = find(theta<0);
theta(ind)=theta(ind)+2*pi;


hold off

%plot(t, theta/pi, '-');

%plot(t/T, amplitude/A0, '-r');

plot(theta/pi, amplitude/A0, '-r');



%----------------------------------------------------
  t = 0:0.0001:tmax;

  %-- A = 20

  %a0 = A0*0.06;  % k01  (good)
  %a0 = A0*0.02;  % k10   (bad)
  %a0 = A0*0.045;   % k11  (good)
  %a0 = A0*0.04;   % k-11 (bad)

  %a0 = A0*0.022;    % k33   (good)
  %a0 = A0*0.013;    % k-33  (wide, but messy)
  %a0 = A0*0.035;    % k03  (good)
  %a0 = A0*0.004;   % k30  (bad)

  %a0 = A0*[0.025, 0.1, 0.2, 0.26, 0.27, 0.4];  % initial amplitudes a0=sqrt(n)
  %a0 = A0*[0.01, 0.025, 0.1, 0.2, 0.4];  % initial amplitudes a0=sqrt(n)
  %a0 = A0*[0.26, 0.27];  % initial amplitudes a0=sqrt(n)

  %-- A = 60
   
   a0 = A0*0.02;  % k01  (good)
  %a0 = A0*0.015;  % k10   (good)
  %a0 = A0*0.01;   % k11  (good)
  %a0 = A0*0.008;  % k-11 (interesting)

  %a0 = A0*0.0006;    % k33   (bad)
  %a0 = A0*0.0008;    % k-33  (bad)
  %a0 = A0*0.006;    % k03  (good)
  %a0 = A0*0.005;   % k30  (good)


  %a0=A0*0.02;
  %a0 = A0*[0.025, 0.1, 0.2, 0.26, 0.27, 0.4];  % initial amplitudes a0=sqrt(n)
  %a0 = A0*[0.025, 0.1, 0.05];  % initial amplitudes a0=sqrt(n)
  %a0 = [1, 2, 3];  % initial amplitudes a0=sqrt(n)



  hold on

  for i = 1:length(a0)
 
     %if i==1 hold off;   else  hold on;  end

     f0 = [pi, (a0(i)).^2];

     f = lsode('fdot', f0, t);

     x = f(:,1); y = f(:,2);

     plot(x/pi,sqrt(y)/A0, 'b'); axis([0,2]);

     %plot(t/T,sqrt(y)/A0, 'b')
     %plot(t,x/pi)

  end

return

  %theta=x; amplitude=sqrt(y); fname="model";

  fid=fopen([fname, '.', fnameext, '.txt'], 'wt');
  fprintf(fid, '# 1.t  2.amplitude/A0  3.theta/pi');
  for i=1:length(tdata)
      fprintf(fid, '%16.6e  %16.6e  %16.6e\n', tdata(i), theta(i)/pi, amplitude(i)/A0);
  end
  fclose(fid);


end

%---------------------

function f = fdot(x,t)

  %A=19.75;
  A=59.73;
  kk=9;

  N = A^2;

  f = [0 0]'; %'

  theta=x(1);
  n=x(2);
  
  f(1) = 2*kk  + 2*(N - 3*n)  + 2*(N - 4*n)*cos(theta);
  f(2) = 2*n*(N - 2*n) * sin(theta);
 
end

%---------------------


