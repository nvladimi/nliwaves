
function psi_images_hpg
%
%   psi_images_hpg
%
%   Make images of Psi magnutude, phase, 
%   and magnitude of gradient of phase.
%   Edit to customize.
%

  for n=1    %:16:257

    one_image(n)

  end

end

%---------------------------------------------------------

function one_image(n)

%-- parameters --

  fbase = 'eps1e-2_4096.psi.'    % name base for input files

  N = 4096;                      % domain is NxN points
  scale=8;                       % reduce images 8 times

%-- filenames -- 

  fname_psi    = [fbase, num2str(n,'%04d')];

  fname_psilin = ['IMG/psilin.',  num2str(n,'%04d'), '.png'];
  fname_psilev = ['IMG/psilev.',  num2str(n,'%04d'), '.png'];
  fname_gphase = ['IMG/gphase.',  num2str(n,'%04d'), '.png'];
  fname_phase =  ['IMG/phase.',   num2str(n,'%04d'), '.png'];

%-- colormaps --

  d=colormap(gray(256));

  b=[(0:63)/63, 1-(0:63)/63, zeros(1,128)];
  g=circshift(b,[1,64]);
  r=circshift(g,[1,64]);
  c=[r',g',b'];
 
  c(1,:) =[0.0, 0.0, 0.0]; 
 
  c( 65,:)=[0.0, 0.0, 0.7]; 
  c(129,:)=[0.0, 0.6, 0.5];
  c(192,:)=[1.0, 0.6, 0.4];
  c(256,:)=[1.0, 1.0, 1.0];

  %r=[zeros(1,128), (1:64)/64, ones(1,64)];
  %g=[zeros(1,64), 0.6*(1:64)/64, 0.6*ones(1,64), 0.4*ones(1,64)+0.6];
  %b=[0.7*(1:64)/64,  -0.2*(1:64)/64+0.7, -0.1*(1:64)/64+0.5, 0.6*(1:64)/64+0.4];

  %q=[(0:127)/128, 1-(0:127)/128];
  %c=[q',q',q'];


%-- read data, compute gradient --

  f = read_psi(fname_psi, N);
  p = angle(f);
  g = phase_grad(p);


%-- rescale data for smaller image

  ind=scale:scale:N;

  p0 =       p(ind, ind) ; 
  g0 = -log( g(ind, ind) );
  f0 =  abs( f(ind, ind) );


%-- amplitude:  continuous grayscale colormap

  sc = imsc(h0,[-1,5]);
  imwrite(sc, d, fname_psilin);

%-- amplitude:  five-color colormap

  sc = imsc(floor(h0), [0,4]);
  imwrite(sc, c, fname_psilev);

%-- phase:  continuous grayscale colormap

  sc  = imsc(p0, [-pi,pi]);
  imwrite(sc, d, fname_phase);

%-- phase gradient:  continuous grayscale colormap

  sc  = imsc(g0, [0,20]);
  imwrite(sc, d, fname_gphase);

end

%---------------------------------------------------------

function g = phase_grad(p)
% p  - phase
% g  - magnitude of gradient of phase

  px = circshift(p,[0,-1])- p;
  py = circshift(p,[-1,0])- p;

  ind1=find(px > pi);
  ind2=find(px <= -pi);

  px(ind1) =  2*pi - px(ind1);
  px(ind2) = -2*pi - px(ind2);
  
  ind1=find(py > pi);
  ind2=find(py <= -pi);

  py(ind1) =  2*pi - py(ind1);
  py(ind2) = -2*pi - py(ind2);

  g = sqrt(px*px + py*py);

end


%---------------------------------------------------------

function sc = imsc(a, limits)

  amin = limits(1);
  amax = limits(2);

  sc   =  uint8(255*(a-amin)/(amax-amin)) + 1;
 
end

%---------------------------------------------------------
