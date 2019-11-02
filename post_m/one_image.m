function one_image(a, amin, amax, fname)
%
%   one_image(a, amin, amax, fname)
%

  a = max(a,amin);
  a = min(a,amax); 

  map = colormap(jet(256));
  %map = colormap(gray(256));

  sc   =  uint8(255*(a-amin)/(amax-amin)) + 1;
  imwrite(flipud(sc), map, fname);
  %saveimage([fname, '.ppm'], flipud(sc), 'ppm', map);

  %imagesc(flipud(a));

end
