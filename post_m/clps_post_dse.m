function clps_post_dse(fname_in, fname_out)
%
%   clps_post(fname_in, fname_out)
%
%   Computes time-dependent collapse data: rescaled time and height,
%   L, V0, beta, gamma, and their derivatives.  Example of usage:
%
%   clps_post('eps1e-2_4096.clps.01' , 'eps1e-2_4096.post.01')


  %R0 = 2.20620086465074607478363406;

  %R0 = 2.2055;       mu =  0;
  %R0 = 2.0333;       mu = -0.3;
  %R0 = 1.7379;       mu = -1.0;
   R0 = 1.7323;       mu = -1.0;

  %-- read data in --

  data = importdata(fname_in);

    t   = data(:, 1);
    h   = data(:, 2);
    phi = data(:, 3);
    ptx = data(:, 4);
    hxx = data(:, 5);
    hyy = data(:, 6);
    pxx = data(:, 7);
    pyy = data(:, 8);

    %-- hack for too many points
    %ind=rem(imax,10):10:length(h);
    %t=t(ind); h=h(ind); ddh=ddh(ind); ddp=ddp(ind);


    %-- compute secondary diagnostics

    L  =  sqrt(h) ./ sqrt(h.^3 + hxx + hyy - mu*h.*ptx );

    L0 = R0./h;
    V0 = h.*L;

    beta  = betaV(V0);
%   gamma = 2*L.^2.*ddp;

    dL   = gradient(L,   t);
    dL0  = gradient(L0,  t);

%   ddL  = gradient(dL,  t);
%   dddL = gradient(ddL, t);

    dB   = gradient(beta, t);

%   ddB  = gradient(dB,   t);
%   dG   = gradient(gamma, t);
%   ddG  = gradient(dG,    t);
    dphi = gradient(phi,  t);

%    beta_tau = dB .* L .* L;


%-- filtering out bad data points --

  ind = find ( abs(dB) < 100 );
  t = t(ind); h=h(ind);  L0=L0(ind); L=L(ind); dL=dL(ind); 
  beta=beta(ind); dB=dB(ind); dphi=dphi(ind); pxx=pxx(ind); pyy=pyy(ind);

  ind = find ( abs(beta) < 1);
  t = t(ind); h=h(ind);  L0=L0(ind); L=L(ind); dL=dL(ind); 
  beta=beta(ind); dB=dB(ind); dphi=dphi(ind); pxx=pxx(ind); pyy=pyy(ind);



%-- output time dependent quantities --

  fid=fopen(fname_out, 'wt');

  fprintf(fid, ... 
    '%% created by \"clps_post_dse.m\" from \"%s\"\n\n', fname_in);
  fprintf(fid, '%%1.t  2.h   3.L0  4.L  5.L_t  6.L0_t  7.beta  8.beta_t  9.phi_t  10.phi_xx  11.phi_yy \n');

  for i=1:length(t) 
      fprintf(fid, ' %17.10e  %17.10e   %13.6e  %13.6e  %13.6e  %13.6e  %13.6e  %13.6e  %13.6e  %13.6e  %13.6e\n',
	      t(i),   h(i),    L0(i),  L(i),  dL(i),  dL0(i),  beta(i),  dB(i),  dphi(i),  pxx(i),  pyy(i));
  end

  fclose(fid);


end

%-----------------------------------------------------


function b = betaV(v)

%  b0   =  4092.03;
%  b1   = -9655.52;
%  b2   =  8540.26;
%  b3   = -3355.19;
%  b4   =   493.905;

% b0  = -1915.04602011082;        #%  +/- 812.7        (42.44%)
% b1  =  4612.09263215859;        #%  +/- 1943         (42.14%)
% b2  = -4165.81973119302;        #%  +/- 1742         (41.82%)
% b3  =  1673.20793610542;        #%  +/- 694.1        (41.49%)
% b4  = -252.235195052472;        #%  +/- 103.7        (41.11%)

% b0  =  149.418061647773           #%  +/- 37.83        (25.32%)
% b1  = -317.84706381942            #%  +/- 90.44        (28.45%)
% b2  =  247.753248058283           #%  +/- 81.08        (32.73%)
% b3  =  -82.4568621578487          #%  +/- 32.3         (39.17%)
% b4  =    9.58977717508506         #%  +/- 4.825        (50.31%)


  b0 =   74.248947577037;        #%  +/- 0.6859       (0.9237%)
  b1 = -138.100177991974;        #%  +/- 1.23         (0.8903%)
  b2 =   86.6016095516639;       #%  +/- 0.7346       (0.8482%)
  b3 =  -18.2552692953927;       #%  +/- 0.1463       (0.8011%)
  b4 =    0;


  b = b0 + b1*v + b2*v.^2 + b3*v.^3 + b4*v.^4;

end

%-----------------------------------------------------
