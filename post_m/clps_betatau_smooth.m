function clps_betatau_smooth

%-------------------
%  inputs
%-------------------

dtau  = 0.1;     % frequency of inputs and width of the averaging window
npoly = 2;       % order of interpolating polynomial


% fname_in='temporal_m277.dat';  fname_out='smooth_m277.dat';
% fname_in='temporal_m2775.dat'; fname_out='smooth_m2775.dat';
% fname_in='temporal_m278.dat';  fname_out='smooth_m278.dat';
% fname_in='temporal_m279.dat';  fname_out='smooth_m279.dat';
% fname_in='temporal_m2795.dat'; fname_out='smooth_m2795.dat';
% fname_in='temporal_m280.dat';  fname_out='smooth_m280.dat';
% fname_in='temporal_m281.dat';  fname_out='smooth_m281.dat';
  fname_in='temporal_m282.dat';  fname_out='smooth_m282.dat';
% fname_in='temporal_m284.dat';  fname_out='smooth_m284.dat';
% fname_in='temporal_m290.dat';  fname_out='smooth_m290.dat';
% fname_in='temporal_m300.dat';  fname_out='smooth_m300.dat';
% fname_in='temporal_m320.dat';  fname_out='smooth_m320.dat';



%--------------------

  %-- read data in --

  data = importdata(fname_in);

  t     = data(:,1);
  tau   = data(:,6);
  L     = data(:,4);
  beta  = data(:,8);

  tmax  = t(end);
  nmax  = fix(tau(end)/dtau) - 1;

  %-- output --

  fid=fopen(fname_out, 'wt');

  fprintf(fid, '# created by \"clps_betatau_smooth.m\" from \"%s\"\n\n', fname_in);
  fprintf(fid, '#1.t  2.tau  3.tmax-t  4.L  5.beta  6.beta_tau\n\n');


  for n=1:nmax

      tau0 = dtau*n;

      ind = find ( tau > dtau*(n-0.5) & tau < dtau*(n+0.5) );

      p = polyfit( tau(ind), beta(ind), npoly );
      beta0 = polyval(p, tau0);
      betatau0 = polyval(polyderiv(p),tau0);

      p = polyfit( tau(ind), L(ind), npoly );
      L0 = polyval(p, tau0);

      p = polyfit( tau(ind), t(ind), npoly );
      t0 = polyval(p, tau0);

      fprintf(fid, ' %22.15e %22.15e %22.15e %22.15e %22.15e %22.15e  %d\n', ...
	      t0, tau0, tmax-t0, L0, beta0, betatau0);
 
   end

  fclose(fid);

end
