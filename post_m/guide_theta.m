N=256; NN=N*N; ink = N/4+1:3*N/4; f=read_psi('n3600a.psi.0000', N);

f=fft2(f)/NN;  a0=arg(f(1,1));  f=fftshift(f);  f=f(ink,ink); f=f(2:end, 2:end); 
n=(abs(f)).^2; a1=arg(f); a2=flipud(fliplr(a1)); th=2*a0 - a1 - a2;

ind=find(th<0); th(ind)=th(ind)+2*pi;
ind=find(th<0); th(ind)=th(ind)+2*pi;

imagesc(th, [0,2*pi])

vtk_any(n, 'nk', 'nk.vtk');
vtk_any(th, 'theta', 'theta.vtk');

