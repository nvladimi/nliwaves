
function vtk_psi
%
%   vtk_psi
%
%   Field data converter.
%   Reduces the size of data and saves it in VTK format for paraviewing. 

  fbaseN = 'eps1e-2_4096.Nblur0000';
  fbaseH = 'eps1e-2_4096.Hblur0000';

  fbase  = 'eps1e-2_4096';

  N = 4096;

  for n=1    %:16:257

    process_one(fbaseN, fbaseH, fbase, N, n)

  end

end

%---------------------------------------------------------

function process_one(fbaseN, fbaseH, fbase, N, n)

%-- rescaling to smaller size --

  scale=4;
  ind=scale:scale:N;
  M = N/scale;



%-- filenames -- 

  fnameN = [fbaseN, num2str(n,'.%04d')];
  fnameH = [fbaseH, num2str(n,'.%04d')];
  fvtkN  = ['VTK/', fbase, '.psi.vtk'];
  fvtkH  = ['VTK/', fbase, '.ener.vtk'];



%-------------- psi ---------------


  fid = fopen(fnameN, 'rb');
  f = fread(fid, N*N, 'float');
  fclose(fid);

  f = reshape(f,N,N);
  f = f(ind, ind); 
  f = sqrt(f);
  f = f';

%-- vtk header --

  fid=fopen(fvtkN, 'wt');
  fprintf(fid, '# vtk DataFile Version 3.0\n');
  fprintf(fid, 'vtk output\n');
  fprintf(fid, 'BINARY\n');
  fprintf(fid, 'DATASET STRUCTURED_POINTS\n');
  fprintf(fid, 'DIMENSIONS %d %d %d\n', M, M, 1);
  fprintf(fid, 'SPACING 1 1 0\n');
  fprintf(fid, 'ORIGIN 0 0 0\n');
  fprintf(fid, 'POINT_DATA %d\n', M*M);
  fprintf(fid, 'SCALARS psi float\n');
  fprintf(fid, 'LOOKUP_TABLE default\n');
  fclose(fid);

%-- vtk data --

  fid=fopen(fvtkN, 'ab', 'ieee-be');
  fwrite(fid, f, 'float');
  fclose(fid);


%-------------- energy ---------------

  fid = fopen(fnameH, 'rb');
  f = fread(fid, N*N, 'float');
  fclose(fid);

  f = reshape(f,N,N);
  f = f(ind, ind); 
  f = f'; 

%-- vtk header --

  fid=fopen(fvtkH, 'wt');
  fprintf(fid, '# vtk DataFile Version 3.0\n');
  fprintf(fid, 'vtk output\n');
  fprintf(fid, 'BINARY\n');
  fprintf(fid, 'DATASET STRUCTURED_POINTS\n');
  fprintf(fid, 'DIMENSIONS %d %d %d\n', M, M, 1);
  fprintf(fid, 'SPACING 1 1 0\n');
  fprintf(fid, 'ORIGIN 0 0 0\n');
  fprintf(fid, 'POINT_DATA %d\n', M*M);
  fprintf(fid, 'SCALARS ener float\n');
  fprintf(fid, 'LOOKUP_TABLE default\n');
  fclose(fid);

%-- vtk data --

  fid=fopen(fvtkH, 'ab', 'ieee-be');
  fwrite(fid, f, 'float');
  fclose(fid);

end

