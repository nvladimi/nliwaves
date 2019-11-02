function clps_history_select


%-- depenpendence on ddh in the long run


  %ih=2;
  %fname_in = 'eps1e-2_4096.clps.16'
  %fname_out = 'clps_hmax-ddh.txt'

%-- depenpendence on ddh in the long run

  ih=5;
  hmax_min = 4;
  fname_in = 'eps1e-2_4096_zoom.bumps.txt'
  fname_out = 'bumps.history.hmax4.txt'

  process_varible(fname_in, fname_out, ih, hmax_min);


end

%-----------------------------------------------------------------

function clps_history_select(fname_in, fname_out, iv, vmin)

  %-- read data in --

  data = importdata(fname_in);

  fid = fopen(fname_out, 'at');

  %-- find begin and end time for each collapse --

  v = data(:,iv);
  ii = find(isnan(v)==1);

  i1 = ii+1;
  i2 = circshift(ii-1, -1); 
  i2(end)=length(v);

  %-- process collapses one-by-one  --

  for n=1:length(i1)

    v   = data(i1(n):i2(n), iv);

    printf('collapse %d\n', n-1);

    imax = find(v==max(v));

    if (imax == length(v))
      printf('  ... abrubted\n');
      continue;
    end

    if (v(imax) < vmax)
      printf('  ... vmax = %f < %f\n', v(imax), vmax);
      continue;
    end

    fprintf(fid, '%12.4e  %12.4e', tmax, hmax);
    fprintf(fid, '%12.4e', vout);
    fprintf(fid, '\n');

  end

  fclose(fid);

end

#-----------------------------------------------------
