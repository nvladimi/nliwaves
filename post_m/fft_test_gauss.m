clear all

%--- input: interval and number of discretization points --

N = 100;

xmin = -10;      
xmax =  10;


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

y = exp( - x.^2 ); 


%-- compute FFT and inverse FFT --

F  = fft(y) ./ (M*N*dk); 

Fi = ifft( F.*M*N*dk );


%-- analytical solution --

f = sqrt(pi)/(2*pi)*exp(-(k.*k)/4);  


%-- plot results --

figure(1); plot(x,Fi,'ro',  x,y, 'b-' );

figure(2); plot(k,F, 'ro', k, f, 'b-');


#--------------

