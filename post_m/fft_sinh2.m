clear all

%--- input: interval and number of discretization points --

N = 4000;

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

y = (sinh(x/2)).^2 ./ ( 1 + sinh(x) ./ x) - abs(x/2); 
ind = find(x==0);
y(ind) = 0;

%-- compute FFT and inverse FFT --

F  = fft(y) ./ (M*N*dk); 

Fi = real(ifft( F.*M*N*dk ));


%-- analytical solution --

%f = sqrt(pi)/(2*pi)*exp(-(k.*k)/4);  


%-- plot results --

%figure(1); plot(x,Fi,'ro',  x,y, 'b-' ); title('y(x)');

%figure(2); plot(k,real(F), 'ro-', k, imag(F), 'bo-'); title('F(k), real=red, imag=blue')

%figure(3); plot(k,real(F), 'ro-'); axis([0.5,0.6, -0.2, 0.2]); grid('on'); title('real F(k), zoom')

figure(4); loglog(k, abs(F),  'bo-',  k, k.^(-2), 'r-');
 axis([1, 100, 1e-4,1]); grid('on'); title('F(k), logscale')

%--------------

