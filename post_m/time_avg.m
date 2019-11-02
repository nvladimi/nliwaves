function time_avg(fname_in, t_index, d_index, t_range)
%
%   function time_avg(fname_in, t_index, d_index, t_range)
%
%   average data in text file,
%   prints average, deviation, and data size on the screen
%
%   Input:
%       fname_in    name base for input files
%       t_index     column containing time (zero - none)
%       d_index     column containing data
%       t_range     time interval (tmin, tmax]  


  d  = importdata(fname_in);

  %dx = 25.6/2048;
  dx = 1;

  if (t_index == 0) 
    i1 = t_range(1)+1;
    i2 = t_range(2);
    d  = d(i1:i2, d_index);
  else
    d = [d(:,t_index), d(:,d_index)]; 

    d  = remove_doubles(d);

    i1 = min( find( d(:, 1) >  t_range(1) ) );
    i2 = max( find( d(:, 1) <= t_range(2) ) );

    d  = d(i1:i2, 2);

  end

  n  = length(d);

  avg = sum(d)/n;
  d = d - avg;
  sigma =  sqrt ( sum(d.*d)/n );

  printf ('\n %g   +/-   %g   (n=  %d )\n\n', avg*dx, sigma*dx, n );


end

%---------------------------------------------------------

function d = remove_doubles(h2)

  [h,i]   = sort(h2(:,1));
  h2(:,1) = h;
  h2(:,2) = h2(i,2);

  for i=2:length(h2)
    if ( h2(i,1)== h2(i-1,1))
      h2(i-1,2) = NaN;
    end
  end

  h=h2(:,2);
  i=isfinite(h);
  h=h(i);
  t=h2(i,1);

  d=[t,h];

end

%---------------------------------------------------------

