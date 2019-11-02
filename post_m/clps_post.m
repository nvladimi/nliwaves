function clps_post(fname_in, fname_out)
%
%   clps_post(fname_in, fname_out)
%
%   Computes time-dependent collapse data: rescaled time and height,
%   L, V0, beta, gamma, and their derivatives.  Example of usage:
%
%   clps_post('eps1e-2_4096.clps.01' , 'eps1e-2_4096.post.01')


  R0=2.20620086465074607478363406;

  %-- read data in --

  data = importdata(fname_in);

  %-- prepare output --

  fid=fopen(fname_out, 'wt');

  fprintf(fid, ... 
    '%% created by \"clps_post.m\" from \"%s\"\n\n', fname_in);
  fprintf(fid, ... 
    '%%1.t  2.h  3.t-tmax  4.(t-tmax)*hmax^2  5.h/hmax  6.L  7.V0  8.gamma  9.beta\n');
  fprintf(fid, ... 
    '%%10.L_t  11.L_tt  12.L_ttt  13.beta_t  14.beta_tt  15.gamma_t  16.gamma_tt\n');
  fprintf(fid, ... 
    '%%17.L0=R0/h  18.L0_t  19.L0_tt  20.dphi\n');


  nout = 0;

  %-- find begin and end time for each collapse --

  tt = data(:,1);
  ii = find(isnan(tt)==1);

  i1 = ii+1;
  i2 = circshift(ii-1, -1); 
  i2(end)=length(tt);

  %-- process collapses one-by-one  --

  for n=1:length(i1)

    fprintf('collapse %d\n', n-1);

    %-- extract data for current collapse --

    t   = tt(i1(n):i2(n));
    h   = data(i1(n):i2(n), 2);
    ddh = data(i1(n):i2(n), 3);
    ddp = data(i1(n):i2(n), 4);
    phi = data(i1(n):i2(n), 5);

    %-- find maximum --

    %if is_duplicate_entry(t)>0 
    %   fprintf('  ... dublicate times\n');
    %   continue
    %end

    imax = find(h==max(h));
    tmax = t(imax);
    hmax = h(imax);
    if (imax == length(h))
      fprintf('  ... abrubted.  Continue post-processing.\n');
      %fprintf(' ... abrubted.  Working on next.\n'); continue;
    end

    %-- hack for too many points
    %ind=rem(imax,10):10:length(h);
    %t=t(ind); h=h(ind); ddh=ddh(ind); ddp=ddp(ind);


    fprintf('  ... hmax = %f\n', hmax);

    %-- compute secondary diagnostics

    V0  = 1./sqrt(1 + 2*ddh./h.^3);
    L   = V0./h;

    beta  = betaV(V0);
    gamma = 2*L.^2.*ddp;

    dL   = gradient(L,   t);
    ddL  = gradient(dL,  t);
    dddL = gradient(ddL, t);

    dB   = gradient(beta, t);
    ddB  = gradient(dB,   t);

    dG   = gradient(gamma, t);
    ddG  = gradient(dG,    t);

    L0   = R0./h;
    dL0  = gradient(L0,   t);
    ddL0 = gradient(dL0,  t);

    dphi = gradient(phi,  t);

    %-- output time dependent quantities

    fprintf(fid, '\n\n%%collapse %d(%d):  t_max= %12.6e    h_max= %12.6e\n', ...
	    nout, n-1, tmax, hmax);
    fprintf(fid, ...
        ' NaN NaN NaN NaN  NaN NaN NaN NaN  NaN NaN NaN NaN  NaN NaN NaN NaN   NaN NaN NaN NaN\n\n');

    for i=1:length(t) 
      fprintf(fid, ...
	' %17.10e %17.10e %13.6e %13.6e %13.6e %13.6e %13.6e',...
	      t(i), h(i), t(i)-tmax, (t(i)-tmax)*hmax*hmax, h(i)/hmax, L(i), V0(i));
      fprintf(fid, ...
	' %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e',...
	      gamma(i), beta(i), dL(i), ddL(i), dddL(i), dB(i), ddB(i), dG(i), ddG(i));
      fprintf(fid, ...
	' %13.6e %13.6e %13.6e %13.6e\n',...
	      L0(i), dL0(i), ddL0(i), dphi(i));
    end

    nout = nout+1;


  end

  fclose(fid);


end

%-----------------------------------------------------


function b = betaVold(v)

  b6=      98.9922359479098;
  b5=   -1403.7339662171987;
  b4=    8286.5470292400041;
  b3=  -26067.2261230988006;
  b2=   46084.0997339975365;
  b1=  -43410.0271983773928;
  b0=   17021.0686019477034;

  b = b0 + b1*v + b2*v.^2 + b3*v.^3 + b4*v.^4 + b5*v.^5 + b6*v.^6;

end


function b = betaV(v)

  b6 =   -5.71952248774065;
  b5 =   64.58276074952414;
  b4 = -286.42012849730742;
  b3 =  610.28639150972106;
  b2 = -581.01029967009026;
  b1 =   96.56198202862554;
  b0 =  131.01891277343202;

  b = b0 + b1*v + b2*v.^2 + b3*v.^3 + b4*v.^4 + b5*v.^5 + b6*v.^6;

end

%-----------------------------------------------------

function imin = ijump(t)

  dtcoarse = t(2)-t(1);
  dt = circshift(t,-1) - t;   
  dt(end)=dt(end-1);

  imin = min( find(dt<0.75*dtcoarse) );

end


%-----------------------------------------------------
