function psi_all
%
%   psi_all
%   Edit script to customize.
%

  fbase='gauss';
  N = 512;
  hmax=10;

  for k=0:4
 
    fname = [fbase, '.psi.', num2str(k,'%04d')];

    psi_images(fname, N, hmax);

  end

end

