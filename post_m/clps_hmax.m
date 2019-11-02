function clps_hmax(fname_in, fname_out, it, ih, iv, h_follow)
%
%   clps_hmax(fname_in, fname_out, it, ih, iv, h_follow)
%
%   Finds "hmax" as function of variable "v" at earlier heights,
%   specified in input array "h_follow".  Examples of usage:
%
%    h_follow=12:2:40;    % values of h, where to measure v
%    it=1;  ih=2;  iv=9;  % indexes of columns, iv=9=beta in *.post 
%    fname_in  = 'eps1e-2_4096.post.01'
%    fname_out = 'eps1e-2_4096.hmax-beta.01'
%
%    h_follow=12:2:40;    % values of h, where to measure v
%    it=1;  ih=2;  iv=3;  % indexes of columns, iv=3=ddh  in *.clps
%    fname_in  = 'eps1e-2_4096.clps.01'
%    fname_out = 'eps1e-2_4096.hmax-ddh.01'


  %-- read data in --

  data = importdata(fname_in);
  
  %-- prepare output --

  fid=fopen(fname_out, 'wt');
  fprintf(fid, '%% created by \"clps_hmax.m\" from \"%s\"\n', fname_in);
  fprintf(fid, '%% input: it = %d, ih = %d, iv = %d, h =', it, ih, iv);
  fprintf(fid, ' %d', h_follow);
  fprintf(fid, '\n\n');

  fprintf(fid, '%%1.tmax  2.hmax  3+.v(h)\n\n');

  %-- find begin and end time for each collapse --

  tt = data(:,it);
  ii = find(isnan(tt)==1);

  i1 = ii+1;
  i2 = circshift(ii-1, -1); 
  i2(end)=length(tt);

  %-- process collapses one-by-one  --

  for n=1:length(i1)

    fprintf('collapse %d\n', n-1);

    %-- extract data for current collapse --

    t   = tt(i1(n):i2(n));
    h   = data(i1(n):i2(n), ih);
    v   = data(i1(n):i2(n), iv);


    %-- find maximum --

    %if is_duplicate_entry(t)>0 
    %   fprintf('  ... dublicate times\n');
    %   continue
    %end

    imax = find(h==max(h));
    tmax = t(imax);
    hmax = h(imax);
    if (imax == length(h))
      fprintf('  ... abrubted\n');
      continue;
    end

    if (hmax < h_follow(1))
      fprintf('  ... hmax = %f < %f\n', hmax, h_follow(1));
      continue;
    end

    fprintf('  ... hmax = %f\n', hmax);


    %-- find v(h0) --

    h1 = min(find( (h_follow<hmax) & h_follow>h(1)));
    h2 = max(find( (h_follow<hmax) & h_follow>h(1)));

    h0 = h_follow(h1:h2);

    t0 = intcubic(h(1:imax), t(1:imax), h0);
    v0 = intcubic(t(1:imax), v(1:imax), t0);

    vout = nan*zeros(1,length(h_follow));
    vout(h1:h2)=v0;

    fprintf(fid, '%12.4e  %12.4e', tmax, hmax);
    fprintf(fid, '%12.4e', vout);
    fprintf(fid, '\n');

  end

  fclose(fid);

end


%-----------------------------------------------------

function y0 = intcubic(x,y, x0)

  coeff = splinecoeff(x,y);

  y0 = zeros(1,length(x0));

  for k=1:length(x0)

    ii=max( find(x<x0(k)) );
    %ii=min( find(x>x0(k)) )-1;

    dx = x0(k)-x(ii);

    %fprintf(' %4d  %16.10f   %16.10f   %16.10f\n', ii, x0(k), x(ii), dx);

    q = coeff(ii,3)*dx;
    q = (q+coeff(ii,2))*dx;
    q = (q+coeff(ii,1))*dx + y(ii);

    y0(k) = q;

  end

end

%-----------------------------------------------------

function coeff=splinecoeff(x,y)
% Calculates coefficients of cubic spline
% Input: x,y vectors of data points
% Output: matrix of coefficients b1,c1,d1;b2,c2,d2;...

  n  = length(x); 
  A  = zeros(n,n);         % matrix A is nxn
  r  = zeros(n,1);

  for i=1:n-1              % define the deltas
    dx(i) = x(i+1) - x(i); 
    dy(i) = y(i+1) - y(i);
  end

  for i=2:n-1              % load the A matrix

    A(i,i-1:i+1) = [dx(i-1), 2*(dx(i-1)+dx(i)), dx(i)];

    r(i) = 3*( dy(i)/dx(i) - dy(i-1)/dx(i-1) ); % right-hand side

  end

  A(1,1) = 1;             % natural spline conditions
  A(n,n) = 1;

% solve for c coefficients
  coeff = zeros(n,3);
  coeff(:,2) = A\r;
              
% solve for b and d
  for i=1:n-1
    coeff(i,3) = (coeff(i+1,2) - coeff(i,2)) / (3*dx(i));
    coeff(i,1) = dy(i)/dx(i) - dx(i)*(2*coeff(i,2) + coeff(i+1,2))/3;
  end

  coeff=coeff(1:n-1,1:3);

end

%-----------------------------------------------------
