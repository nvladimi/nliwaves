
function psi_images_NH
%
%   psi_images)NH
%
%   Make images of blurred energy and number of particles for
%   movies with circles. Circles are superimposed by python script.
%   Edit to customize.
%

  for n=1    %:16:257

    one_image(n);

  end

end

%---------------------------------------------------------

function one_image(n)
%
% n - file number and image number


%-- parameters --

  fbase = 'eps1e-2_4096.psi.';   % base name for input files

  N = 4096;                      % domain is NxN points
  L = 25.6;                      % length of the domain
  scale=8;                       % reduce images 8 times

  rr1 = 0.3^2;                   % blurring radii
  rr2 = 0.6^2;
  rr3 = 1.2^2;
  rr4 = 2.4^2;    % 2.73/1.1;

%-- filenames -- 

  fname_psi   = [fbase, num2str(n,'%04d')];

  fname_E0  = ['IMG/en00.',  num2str(n,'%04d'), '.png'];
  fname_E1  = ['IMG/en03.',  num2str(n,'%04d'), '.png'];
  fname_E2  = ['IMG/en06.',  num2str(n,'%04d'), '.png'];
  fname_E3  = ['IMG/en12.',  num2str(n,'%04d'), '.png'];
  fname_E4  = ['IMG/en24.',  num2str(n,'%04d'), '.png'];

  fname_NP0 = ['IMG/np00.',  num2str(n,'%04d'), '.png'];
  fname_NP1 = ['IMG/np03.',  num2str(n,'%04d'), '.png'];
  fname_NP2 = ['IMG/np06.',  num2str(n,'%04d'), '.png'];
  fname_NP3 = ['IMG/np12.',  num2str(n,'%04d'), '.png'];
  fname_NP4 = ['IMG/np24.',  num2str(n,'%04d'), '.png'];

%-- colormaps --

  d=colormap(gray(256));

  b=[(0:63)/64, 1-(0:63)/64, zeros(1,128)];
  g=circshift(b,[1,64]);
  r=circshift(g,[1,64]);
  c=[r',g',b'];
 
  c(1,:) =[0.0, 0.0, 0.0]; 
 
  c( 65,:)=[0.0, 0.0, 0.7]; 
  c(129,:)=[0.0, 0.6, 0.5];
  c(192,:)=[1.0, 0.6, 0.4];
  c(256,:)=[1.0, 1.0, 1.0];

%-- read data, compute energy --

  f = read_psi(fname_psi, N);
  u = real(f);
  v = imag(f);

  ux = circshift(u,[0,-1]) - circshift(u,[0,1]);
  uy = circshift(u,[-1,0]) - circshift(u,[1,0]);
  vx = circshift(v,[0,-1]) - circshift(v,[0,1]);
  vy = circshift(v,[-1,0]) - circshift(v,[1,0]);

  dx = L/N;
  over4dxdx= 1/(4*dx*dx);

  E  = (ux.*ux + uy.*uy + vx.*vx + vy.*vy)*over4dxdx;
  NP = u.*u + v.*v;
  E = E - NP.*NP * 0.5;

  clear u; clear v; clear ux; clear uy; clear vx; clear vy;
 
%-- prepare gaussian(s) for convolution

  x=linspace(-L/2, L/2, N+1);
  x=x(1:N);

  [x,y]=meshgrid(x);
  rr = x.*x + y.*y;
  rr = circshift(rr, [N/2,N/2]);

  G = exp ( - rr/rr1 ) / (pi*rr1);  fG1 = fftn(G);
  G = exp ( - rr/rr2 ) / (pi*rr2);  fG2 = fftn(G);
  G = exp ( - rr/rr3 ) / (pi*rr3);  fG3 = fftn(G);
  G = exp ( - rr/rr4 ) / (pi*rr4);  fG4 = fftn(G);

  clear x; clear y; clear rr, clear G;

%-- blur the energy and number of particles

  ind=scale:scale:N;

  NP0 = real(NP(ind, ind));
  fNP = fftn(NP)/(2*pi*N);
  clear NP;

  f = fG1.*fNP;  f = ifftn(f);   NP1 = real(f(ind, ind));
  f = fG2.*fNP;  f = ifftn(f);   NP2 = real(f(ind, ind));
  f = fG3.*fNP;  f = ifftn(f);   NP3 = real(f(ind, ind));
  f = fG4.*fNP;  f = ifftn(f);   NP4 = real(f(ind, ind));


  E0 = real(E(ind, ind));
  fE = fftn(E)/(2*pi*N);
  clear E;

  f = fG1.*fE;  f = ifftn(f);   E1 = real(f(ind, ind));
  f = fG2.*fE;  f = ifftn(f);   E2 = real(f(ind, ind));
  f = fG3.*fE;  f = ifftn(f);   E3 = real(f(ind, ind));
  f = fG4.*fE;  f = ifftn(f);   E4 = real(f(ind, ind));


%-- images: continuous grayscale colormap

  sc = imsc(E0, [-30,70]);
  imwrite(sc, d, fname_E0);

  sc = imsc(E1, [-30,70]);
  imwrite(sc, d, fname_E1);

  sc = imsc(E2, [-30,70]);
  imwrite(sc, d, fname_E2);

  sc = imsc(E3, [-30,70]);
  imwrite(sc, d, fname_E3);

  sc = imsc(E4, [-30,70]);
  imwrite(sc, d, fname_E4);


  sc = imsc(NP0, [-2,10]);
  imwrite(sc, d, fname_NP0);

  sc = imsc(NP1, [-2,10]);
  imwrite(sc, d, fname_NP1);

  sc = imsc(NP2, [-2,10]);
  imwrite(sc, d, fname_NP2);

  sc = imsc(NP3, [-2,10]);
  imwrite(sc, d, fname_NP3);

  sc = imsc(NP4, [-2,10]);
  imwrite(sc, d, fname_NP4);


end

%--------------------------------------------------

function sc = imsc(a, limits)

  amin = limits(1);
  amax = limits(2);

  sc   =  uint8(255*(a-amin)/(amax-amin)) + 1;
 
end

%--------------------------------------------------
