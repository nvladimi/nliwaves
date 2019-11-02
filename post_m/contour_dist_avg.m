function contour_dist_avg(fbase_in, fname_out, tmax, nfiles)

  dt = 1;

  X=[];
  Y=[];

  for n = nfiles

    fname_in = [fbase_in, '.',  num2str(n,'%04d')];

    r = importdata(fname_in);

    s = r - circshift(r,[1,;]);

    s(1,:) = 0;

    s = sqrt( s(:,1).*s(:,1) + s(:,2).*s(:,2) ); 

    for i=2:length(s)
	   s(i) = s(i) + s(i-1);
    end

    t = dt:dt:s(end,1);

    x = interp1(s, r(:,1), t);
    y = interp1(s, r(:,2), t);

    clear s
    clear r

    m = fix( length(t)/tmax );
    
    x = reshape(x(1:m*tmax), tmax, m); 
    y = reshape(y(1:m*tmax), tmax, m); 

    X = [X,x]; 
    Y = [Y,y]; 

  end

  [n,m] = size(X);

  X0 = repmat(X(1,:), n, 1);
  X  = X - X0;

  Y0 = repmat(Y(1,:), n, 1);
  Y  = Y - Y0;

  RR = X.*X + Y.*Y;

  rr = sum(RR,2)/m;


  %-- output -- 

  fid = fopen(fname_out,'wt');

  fprintf(fid, '%% created by \"contour_dist_avg.m\" from \"%s.*\"\n', fbase_in);
  fprintf(fid, '%% 1.s  2.r^2  (averaged over %d segments)\n\n', m);

  for i=1:n
     fprintf(fid,' %16.8e  %16.8e\n', t(i), rr(i) );
  end

  fclose(fid);



end

%---------------------------------------------------------
