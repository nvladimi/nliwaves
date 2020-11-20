function myrand48_seed(fname, np, seed)
%
%  myrand48_seed(fname, np, seed)
%
%  Input:
%       fname:  file name to create, in format "*.seed.????"
%       np:     number of processes to run on
%       seed:   master seed (interger)
%

  jmax = 100;
  
  rand ("twister", seed);

  j = rand(np*3, 1); 
  j = fix(j*jmax)+1;

  ind = zeros(np*3,1);
  ind(1) = j(1);

  for i=2:np*3
      ind(i) = ind(i-1) + j(i);
  end
  
  x = rand(ind(end), 1);
  x = x(ind);

  xmax = 2^16 - 1;
  y = fix(x*xmax); 
  z = uint16(y);

  fid = fopen(fname, 'wb');
  fwrite(fid, z, 'uint16');
  fclose(fid);
  
end

