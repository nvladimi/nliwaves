function f = guide_forcing(forcing_type, N)
% N is number of dealiased modes.

  if (forcing_type == 2)
    f=forcingDF(N);
  end

  if (forcing_type == 3)
    f=forcingHiK(N);
  end

end


%--------------------------------------------

function f=forcingDF(N)

  s = N/128;

  n  = 128*s;
  kl =  28*s;
  kr =  32*s;
  kd =  42*s;

  a = 6.4e-2/(s*s);
  b = 10*a;
 
  %-- forcing shell --
 
  kkl = kl*kl;
  kkr = kr*kr;


  %-- grid --

  k = [-n/2:n/2-1]; 

  [kx, ky] = meshgrid (k, k);

  kk = kx.^2 + ky.^2;


  %-- pumping --

  pump = (kk - kkl).*(kkr - kk);  
  pump = max(pump, 0);
  pump = sqrt(pump);

  %-- damping --

  q  = sqrt(kk)/kd;
  qq = q.*q;


  h  = zeros(n,n);

  ind = find(q>=1);
  h(ind) = 1 - 5/6. * exp(0.5 - 0.5*qq(ind)); 

  ind = find(q<1);
  h(ind) = exp(5 - 5./qq(ind)) ./ (6*qq(ind).*q(ind).*qq(ind));

  ind = find(q<1.e-6);
  h(ind) = 0;

  damp = kk.*h;


  %-- total ---

  f = a*pump - b*damp;

  ind = find(kk<kkl); 
  f(ind) = 0;

end

%--------------------------------------------

function f=forcingHiK(N)

  s = N/128;

  n  = 128*s;
  kl =  28*s;
  kr =  32*s;

  a = 4;
 
  %-- forcing shell --
 
  kkl = kl*kl;
  kkr = kr*kr;


  %-- grid --

  k = [-n/2:n/2-1]; 

  [kx, ky] = meshgrid (k, k);

  kk = kx.^2 + ky.^2;


  %-- forcing --

  q = (sqrt(kk)-kl)/(kr-kl); 

  f =  4*a*q.*(1 - q);

  ind = find(kk<=kkl); 
  f(ind) = 0;

end

%--------------------------------------------
