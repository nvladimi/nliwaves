function psi_contours(fnum)

snum = ['.', num2str(fnum, '%04d')];

fname_in = ['m4_a1e2_0005_f1e2.psi', snum] ;  
fbase1 = ['contours/m4_a1e2_f1e2.contour', snum];
fbase2 = ['contours/m4_a1e2_f1e2.contraw', snum];


%fname_in = 'm2_a1e2_f1e2.psi.0010';  fbase = 'contour_test';

N = 4096;  M = 304;

Lmin = N/4;

%-- read field and extract contours --

  f=read_psi(fname_in, N);

  ind = M+2:N-M+1;
  q=fftn(f); q(ind,:)=0; q(:,ind)=0; q=angle(ifftn(q));

  % ind=2:2:N; q=q(ind,ind); N=N/2;  % testing 

  c = contourc(q, [0,0]);

  fid=fopen([fbase1, '.out'], 'wt');
  fprintf(fid, '%% created by \"psi_contours.m\" from \"%s\"\n', fname_in);
  fprintf(fid, '%% data format:  1.x 2.y \n\n');
  fclose(fid);

%-- select countours for output --

  nout = 0;
  i0   = 1;

  while (i0<length(c)) 

    n = fix(c(2,i0));
    i1 = i0+1;
    i2 = i0+n;
    i0 = i2+1;
     

    c0 = c(:,i1:i2); 
    dc = circshift(c0,[0,-1]) - c0;
    dc = dc(:,1:end-1);
    dc = sqrt(dc(1,:).^2 +  dc(2,:).^2);
    L  =  sum(dc);

    if (L < Lmin) 
        continue;
    end

    %-- uncentered (raw) output --

    fname_out = [fbase2, '.', num2str(nout, '%04d')];
    fid=fopen(fname_out, 'wt');
    for i=1:length(c0)
        fprintf(fid, ' %13.6e %13.6e\n', c0(1,i), c0(2,i));
    end
    fclose(fid);



    %-- recenter at zero ---

    x = c0(1,:);
    y = c0(2,:); 
    x0 = x(1);
    y0 = y(1);

    if (y0 == 1)          % crossing bottom
  
       c0(1,:) = x - x0;
       c0(2,:) = y - y0;

    elseif (x0 == 1)      % crossing left 

       c0(1,:) = y0 - y;
       c0(2,:) = x - x0;

    elseif (y0 == N)      % crossing top

       c0(1,:) = x0 - x;
       c0(2,:) = y0 - y;

    elseif (x0 == N)      % crossing right

       c0(1,:) = y - y0;
       c0(2,:) = x0 - x;

    else                  % loop

       [m,i] = min(y);
       c0(1,:) = x - x(i);
       c0(2,:) = y - y(i);

       c0 = circshift(c0,[0,-i+1]);
   
    end

    %-- output --

    fname_out = [fbase1, '.', num2str(nout, '%04d')];
    fid=fopen(fname_out, 'wt');
    for i=1:length(c0) 
	fprintf(fid, ' %13.6e %13.6e\n', c0(1,i), c0(2,i));
    end
    fclose(fid);

    fid=fopen([fbase1, '.out'], 'at');
    fprintf(fid, '%% contour %d:  n= %d   L=%13.6e \n', nout, n, L);
    fclose(fid);

    nout = nout+1;

  end

end

