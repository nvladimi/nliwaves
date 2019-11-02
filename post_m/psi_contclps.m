function psi_contclps(fbase, fnum)

snum = num2str(fnum, '%04d');

fname_in = [fbase,'.psi.', snum] ;  

fbase1 = ['post/', fbase, '.cAmp.', snum];
fbase2 = ['post/', fbase, '.cH.',   snum];
fbase3 = ['post/', fbase, '.pAmp.', snum];


N = 320; L=32.0; dx = L/N;

AmpVal = [0.25, 0.5, 1.0, 2.0, 4.0, 8.0];

x=-L/2:L/N:L/2-L/N; k=[0:N/2-1  -N/2:-1]*2*pi/L;

[kx,ky]=meshgrid(k,k);

dx=L/N;

Lmin=0;


%-- read field and extract contours --

  psi = read_psi(fname_in, N);

  f = fftn(psi);
  px = ifftn(1i*kx.*f);
  py = ifftn(1i*ky.*f);

  %px1 = (circshift(psi,[0,-1]) - circshift(psi,[0,1]))/(2*dx);
  %py1 = (circshift(psi,[-1,0]) - circshift(psi,[1,0]))/(2*dx);

  f = abs(psi);
  H = px.*conj(px) + py.*conj(py) - 0.25*f.^4;

% figure(1, 'name', 'amplitude'); imagesc(f); axis('off', 'square')
% figure(2, 'name', 'amplitude'); imagesc(H, [-0.1,0.1]); axis('off', 'square')


%-- profiles of Amplitude --

  fmax=max(max(f));

  fx = f(:, N/2+1);
  fy = f(N/2+1, :);

  fid=fopen( fbase3, 'wt');
  fprintf(fid, 
          '%% profiles of amplitude created by \"psi_contclps.m\" from \"%s\"\n', fname_in);
  fprintf(fid, '%% data format:  1.x,y   2.amp(x)   3.amp(y)\n\n');
 
 
  for i=1:N
     fprintf(fid, ' %13.6e  %13.6e  %13.6e\n', (i-N/2-1)*dx, fx(i), fy(i) );
  end
 
  fclose(fid);


%-- contours of amplitude --

%for v = 1:length(AmpVal);
%  val=AmpVal(v);

  val = fmax/2;
  c = contourc(f, [val, val]);

  if length(c) == 0 
      continue;
  end

  %fname_out = [fbase1, '.', num2str(v, '%04d')];
  fname_out = fbase1;

  fid=fopen(fname_out, 'wt'); 
  fprintf(fid, 
    '%% contour |psi| = %f created by \"psi_contclps.m\" from \"%s\"\n', 
     val, fname_in);
  fprintf(fid, '%% data format:  1.x 2.y \n\n');

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

    for i=1:length(c0)
        fprintf(fid, ' %13.6e %13.6e\n', c0(1,i), c0(2,i));
    end

    fprintf(fid, '\n\n');

  end

  fclose(fid);

%end


%-- contours of zero Hamiltonian --

  c = contourc(H, [0,0]);

  fid=fopen( fbase2, 'wt');
  fprintf(fid, '%% countour H = 0 created by \"psi_contclps.m\" from \"%s\"\n', fname_in);
  fprintf(fid, '%% data format:  1.x 2.y \n\n');
 
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

    for i=1:length(c0)
        fprintf(fid, ' %13.6e %13.6e\n', c0(1,i), c0(2,i));
    end
    fprintf(fid, '\n\n');

  end
 
  fclose(fid);

end

