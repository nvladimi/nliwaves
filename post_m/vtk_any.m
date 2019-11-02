
function vtk_any(v, vname, fname)

  [N,M] = size(v);

%-- vtk header --

  fid=fopen(fname, 'wt');
  fprintf(fid, '# vtk DataFile Version 3.0\n');
  fprintf(fid, 'vtk output\n');
  fprintf(fid, 'BINARY\n');
  fprintf(fid, 'DATASET STRUCTURED_POINTS\n');
  fprintf(fid, 'DIMENSIONS %d %d %d\n', N, M, 1);
  fprintf(fid, 'SPACING 1 1 1\n');
  fprintf(fid, 'ORIGIN 0 0 0\n');
  fprintf(fid, 'POINT_DATA %d\n', N*M);
  fprintf(fid, 'SCALARS %s float\n', vname);
  fprintf(fid, 'LOOKUP_TABLE default\n');
  fclose(fid);

%-- vtk data --

  fid=fopen(fname, 'ab', 'ieee-be');
  fwrite(fid, v, 'float');
  fclose(fid);

end

