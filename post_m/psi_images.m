function psi_images(fname, N, hmax)
%
%   psi_images(fname, N, hmax)
%
%   Make images of Psi magnutude and phase. 
%
%   Input:  
%
%       fname       binary input file
%       N           size of simulation
%       hmax        maximum magnitude to scale colormap
%


  map = colormap(jet(256));

%-- filenames -- 

  fname_psi    = [fname, '.h.png'];
  fname_phase  = [fname, '.p.png'];

  f = read_psi(fname, N);

%-- amplitude --

  a    =  abs(f);

  amin = 0;
  amax = hmax; 

  a = max(a,amin);
  a = min(a,amax); 

  sc   =  uint8(255*(a-amin)/(amax-amin)) + 1;
  imwrite(sc, map, fname_psi);

%-- phase --

  a    =  angle(f);

  amin = -pi;
  amax =  pi; 

  a = max(a,amin);
  a = min(a,amax); 

  sc   =  uint8(255*(a-amin)/(amax-amin)) + 1;
  imwrite(sc, map, fname_phase);

end


%---------------------------------------------------------
