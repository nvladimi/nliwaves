function clps_all
%
%   clps_all
%   Edit script to customize.
%


  %base='eps1e-2_4096';
  base = 'seed1b';

  h_follow = 16:4:40;    % values of h, where to measure v

  % indexes of columns

  it    = 1;
  ih    = 2;  
  ibeta = 9;
  iddh  = 3;

 
  for k=1:15
 
    fname_clps = [base, '.clps.',      num2str(k,'%02d')];
    fname_post = [base, '.post.',      num2str(k,'%02d')];
    fname_beta = [base, '.hmax-beta.', num2str(k,'%02d')];
    fname_ddh  = [base, '.hmax-ddh.',  num2str(k,'%02d')];

    clps_post(fname_clps, fname_post)
    clps_hmax(fname_post, fname_beta, it, ih, ibeta, h_follow)
    clps_hmax(fname_clps, fname_ddh,  it, ih, iddh,  h_follow)

  end

end

