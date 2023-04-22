function giaflucs_profiles(fbase)

%-- read pararemters and decide to read seeds from file or scan consecutively --

  [N, dx, nseeds, I1, I2, factor, numtheta] = read_param_profiles(fbase);

  favg =   [fbase, '_ravg.dat'];
  fseeds = [fbase, '_seeds.dat'];

  if (nseeds>0)
     isconsecutive = 1;
     saveseeds=[];
  else
     load(fseeds);
     %seeds = seeds(1:1);  %debugging
     nseeds = length(seeds);
     isconsecutive = 0;
  end
       
%-- grid for interpolation and averaging --- 

  Nfine =   N*factor;    % size of averaged arrays
  dxfine = dx/factor;    % resolution of averaged arrays

%-- fisrt round: computing average 2D profile --

  Iavg = zeros(Nfine, Nfine);

  for s=1:nseeds

     if isconsecutive
        seed = s;
     else
        seed = seeds(s);
     end 

     U = giaflucs(fbase, seed);

     if isconsecutive
        Imax = max(max(U.*conj(U)));
        if (Imax > I1 && Imax <= I2)
	  saveseeds = [saveseeds;seed];
          Iavg = Iavg + reCenter(U, factor);
        end
     else
          Iavg = Iavg + reCenter(U, factor);
     end 
  end

  if isconsecutive
     seeds = saveseeds;
     nseeds = length(seeds);
  end

  Iavg = Iavg/nseeds;


%-- make angle-averaged profile and 2D field --

  Ir = angle_average(Iavg, dxfine, numtheta);

  Iavg = radial_2D(Ir, dxfine, Nfine);


%-- second round: computing 2D sigma --

  sigmasq = zeros(Nfine, Nfine);

  for s=1:nseeds

     seed = seeds(s);
 
     U = giaflucs(fbase, seed);

     diff = reCenter(U, factor) - Iavg;

     sigmasq =  sigmasq + diff.*diff;

  end

  sigmasq = sigmasq/nseeds;


%-- make angle-averaged profile --

  sigma = sqrt (angle_average(sigmasq, dxfine, numtheta) );


%-- save to files --

  if isconsecutive
     save(fseeds, "seeds");
  end

  r = (0:length(Ir)-1)*dxfine;
  fwhm = 2*interp1(Ir, r, Ir(1)*0.5, 'cubic');

  fid=fopen(favg, 'wt');
 
  fprintf(fid, '# intensity profile averaged over angle and over realization\n');
  fprintf(fid, '# Imax =  %12.6e   FWHM =  %12.6e \n', Ir(1), fwhm);
  fprintf(fid, '#_1.r  2.avg(I)   3.sigma(I) \n\n');

  for i=1:length(Ir)
	  fprintf(fid, ' %12.6e   %12.6e   %12.6e\n', r(i), Ir(i), sigma(i) );
  end

end


%------------------------------------------------------

function I = reCenter(U, factor)

    N  = length(U);

    %figure(1); imagesc(U.*conj(U)); axis('equal');   %debuging

    %-- interpalate whole array to finer grid --

    N1 = N;
    N2 = N*factor;

    fp1 = fft2(U);
    fp1 = fftshift(fp1);

    ind = (N2/2 - N1/2 + 1) : (N2/2 + N1/2);

    fp2 = zeros(N2,N2);
    fp2(ind,ind) = fp1;
    fp2 = ifftshift(fp2);
    U = (N2/N1)**2 * ifft2(fp2);

    I  = U.*conj(U);

    %figure(2); imagesc(I); axis('equal');   %debuging

    %-- find maximum of intensity and recenter --

    ind = find(I == max(max(I)));
    [i1, i2] = ind2sub ([N2, N2], ind);
    I = circshift(I, [N2/2-i1, N2/2-i2]);
 
    %figure(3); imagesc(I); axis('equal');  % debuging

    %ind = find(I == max(max(I))); [i1, i2] = ind2sub ([N2, N2], ind) %debugging

end

%------------------------------------------------------

function dat_avg = angle_average(dat, dx, numtheta)
 
  dr = dx; 

  N = length(dat);
  
  x1 = ( -N/2+1:N/2)*dx;
  y1 = x1;
  [x, y] = meshgrid (x1, y1);

  %ind = find(dat == max(max(dat))); xmax=x(ind), ymax=y(ind) %debugging

%-- interpolation --

  thetaI = 2*pi*(0:numtheta-1)/numtheta; 
  rI = (0:N/2-1)*dr;

  [theta,r] = meshgrid(thetaI,rI);

  xI = r.*cos(theta);
  yI = r.*sin(theta);

  dat_cyl = interp2(x,y, dat, xI,yI, 'cubic');

%-- averaging -- 

  dat_avg = sum(dat_cyl,2)/numtheta;

end

%------------------------------------------------------

function dat = radial_2D(dat_rad, dx, N)

  dr = dx; 
 
  r1D = (0:N-1)*dr;
  dat1D = r1D*0;
  dat1D(1:length(dat_rad)) = dat_rad;

  x1 = ( -N/2+1:N/2)*dx;
  y1 = x1;
  [x, y] = meshgrid (x1, y1);
  r2D = sqrt(x.*x + y.*y); 

  dat = interp1(r1D, dat1D, r2D, 'cubic');

  %ind = find(dat == max(max(dat))); xmax=x(ind), ymax=y(ind) %debugging

end

%------------------------------------------------------


function [N, dx, nseeds, I1, I2, factor, numtheta] = read_param_profiles(fbase)

  fid = fopen([fbase, '.param']);
  if strncmp(version,'3.4.0',5)
       s   = textscan(fid,'%s%s', 'commentstyle', 'shell');
  else
       s   = textscan(fid,'%s%s', 'commentstyle', '#');
  end
  fclose(fid);

  %----------------------

  ind = strncmp('N',s{1},1);      
  if (ind<0) disp('undefined "N"');
  else a=s{2}(ind){1}; N =      str2num(a); 
  end

  ind = strncmp('dx',s{1},2);
  if (ind<0) disp('undefined "dx"');  
  else a=s{2}(ind){1}; dx =     str2num(a);
  end

  ind = strncmp('nseeds',s{1},6);  
  if (ind<0) disp('undefined "nseeds"');
  else a=s{2}(ind){1}; nseeds =  str2num(a);
  end

  ind = strncmp('I1',s{1},2);     
  if (ind<0) disp('undefined "I1"');
  else a=s{2}(ind){1}; I1 =     str2num(a);
  end

  ind = strncmp('I2',s{1},2); 
  if (ind<0) disp('undefined "I2"');
  else a=s{2}(ind){1}; I2 = str2num(a);
  end

  ind = strncmp('factor',s{1},6); 
  if (ind<0) disp('undefined "factor"');
  else a=s{2}(ind){1}; factor = str2num(a);
  end

  ind = strncmp('numtheta',s{1},8); 
  if (ind<0) disp('undefined "numtheta"');
  else a=s{2}(ind){1}; numtheta = str2num(a);
  end

end




%------------------------------------------------------
