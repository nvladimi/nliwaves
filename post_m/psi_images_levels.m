function psi_imageslevels


  printf('#1.time  2.hmax  3.level  4.radius  5.havg  6.N\n\n');

    for n=13    %:16:257

    %one_image(n)
    one_analysis(n)

  end

end

%---------------------------------------------------------

  function one_image(n)

    fbase = 'ga5x5_001.psi.';    % name base for input files

  N = 512;                      % domain is NxN points
  scale=2;                      % reduce images 8 times

  %levels = [0.5,1, 2,3,3,3]
  levels = 0.5:0.5:2;

%-- filenames -- 

  fname_psi    = [fbase, num2str(n,'%04d')];
  fname_psilev = ['IMG/psilev.',  num2str(n,'%04d'), '.png'];


%-- read data, compute gradient --

  ind=scale:scale:N;

  f = read_psi(fname_psi, N);
  h =  abs( f(ind, ind) );

  h0=uint8(zeros(size(h))); 

  nlevels=length(levels);

  for n=1:nlevels

      ind = find( h > levels(n));
      h0(ind) = n;

  end

 %imagesc(h0); axis('equal');

  c = colormap(gray( nlevels+1 ));
  imwrite(h0+1, c, fname_psilev);


end

%---------------------------------------------------------
  function one_analysis(n)

  fbase = 'ga5x5_001.psi.';    % name base for input files
  N = 512*2;                     % domain is NxN points
  L = 12.8;
  scale = 2*2;
  
  levels = [1.5,2];
  %levels = 1;
   max_recursion_depth(1024);

%-- filenames -- 

  dx = L/N*scale;

  fname_psi = [fbase, num2str(n,'%04d')];

  ind=scale:scale:N;

  f = read_psi(fname_psi, N);
  h =  abs( f(ind, ind) );

  for lev=1:length(levels)

    H=logical(zeros(size(h)));
    ind = find( h > levels(lev));
    H(ind) = 1;

    S = detect(H);

    imagesc(S); axis('equal'); pause(1.0);

    Smax=max(max(S));

    for iS=1:Smax 
      ind = find(S==iS);
      A = length(ind);
      Ravg = sqrt(A)/pi * dx; 
      havg = sum(sum(h(ind)))/A;
      hmax = max(max(h(ind)));
      Ntot = sum(sum(h(ind).^2))*dx*dx;

      printf('%3d  %f  %f  %f  %f  %f\n',
	     n, hmax, levels(lev), Ravg, havg, Ntot);
    end


  end



end


%---------------------------------------------------------
