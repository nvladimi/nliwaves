function Y=contour_map_avg(fbase_in, fname_out, imax, nfiles)

  dx = 1;
  %imax = 2000;

  Y=[];

  for n = nfiles

    fname_in = [fbase_in, '.',  num2str(n,'%04d')];

    d = importdata(fname_in);

    x = dx:dx:d(end,1);

    y = interp1(d(:,1), d(:,2), x);

    m = fix( length(y)/imax );
    
    y = reshape(y(1:m*imax), imax, m); 

    Y = [Y,y]; 

  end


  [n,m] = size(Y);

  Y0 = repmat(Y(1,:), n, 1);

  Y = Y - Y0;

  D = sum(Y.*Y,2)/m;
  T =  dx*transpose(0:n-1); 


  %-- output -- 

  fid = fopen(fname_out,'wt');

  fprintf(fid, '%% created by \"contour_map_avg.m\" from \"%s.*\"\n', fbase_in);
  fprintf(fid, '%% 1.T  2.D^2  (averaged over %d segments)\n\n', m);

  for i=1:n
     fprintf(fid,' %16.8e  %16.8e\n', T(i), D(i) );
  end

  fclose(fid);


end

%---------------------------------------------------------
