clear all

%--- input: interval and number of discretization points --

N = 400;

xmin = -100;      
xmax =  100;


%---- supplemental arrays:  discretization --

L  = xmax-xmin;
dx = L/N; 
x  = xmin + (0:N-1)*dx; 

%---- supplemental arrays:  integer modes --

m=0:N-1;
m0 = floor((N-1)/2);   
m = circshift(m-m0,[0,-m0]);  
M = (-1).^m;


%---- supplemental arrays:  "physical" modes --

dk = 2*pi/L;
k = dk * m;            


%-- define function ---

y = sinh(x/2) ./ ( 1 + sinh(x) ./ x); 
ind = find(x==0);
y(ind) = 0;

%-- compute FFT and inverse FFT --

F  = fft(y) ./ (M*N*dk); 

Fi = real(ifft( F.*M*N*dk ));


%-- analytical solution --

%f = sqrt(pi)/(2*pi)*exp(-(k.*k)/4);  


%-- plot results --

figure(1); plot(x,Fi,'ro',  x,y, 'b-' );

figure(2); plot(k,real(F), 'ro-', k, imag(F), 'bo-');

figure(4); plot(k,imag(F), 'b-'); axis([1,1.5, -1e-3,1e-3]); grid('on');

figure(3); semilogy(k, abs(F),  'bo-'); axis([0, 5, 1e-10,1]); grid('on');

#--------------

