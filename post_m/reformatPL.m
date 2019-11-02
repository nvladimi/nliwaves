function reformatPL(fbase)
%
%   reformatPL(fbase)
%
%   reads in *.clps and *.dat files, writes out text file 
%   in PL's format *_pl.txt only for isolated collapses
%

  %-- read data in --

  clps = importdata([fbase, '.clps']);
  data = importdata([fbase, '.dat']);

  %-- match two sets of times --

  clps = clps(2:end,:);

  t1 = clps(:,1);
  t2 = data(:,1);
  tbeg = t1(1);
  ibeg = find(t2 == tbeg);
  tend = t1(end);
  iend = max(find(t2 <= tend));

  data = data(ibeg:iend,:);
  t2 = data(:,1);

  index = zeros(length(t2),1);

  for i2=1:length(t2)
    i1 = find(t1 == t2(i2));
    index(i2) = i1;
  end

%-- data from collapse file --

  h   = clps(:,2);
  ddh = clps(:,3);
  ddp = clps(:,4);
  p   = clps(:,5);
  u   = h.*cos(p);
  v   = h.*sin(p);


%-- arrays for N and energies

  N  = zeros(length(t1),1);
  Ek = zeros(length(t1),1);
  Ep = zeros(length(t1),1);

  N(index)  =  data(:,4);
  Ek(index) =  data(:,5);
  Ep(index) = -data(:,6);
 

%-- output
 
  fid=fopen([fbase, '_pl.txt'], 'wt');

  %fprintf(fid, '%%1.t  2.|psi|  3.N   4.Ek  5.Ep  6.empty  7.empty  8.empty  9.phi');
  %fprintf(fid, '  10.re(psi)  11.im(psi)  12.empty  13.|psi|_rr  14.phi_rr\n');

  for i=1:length(t1)

    fprintf(fid, '%14.6e %14.6e %14.6e %14.6e %14.6e  0  0  0 ',...
	    t1(i), h(i), N(i), Ek(i), Ep(i));
    fprintf(fid, ' %14.6e %14.6e %14.6e  0  %14.6e %14.6e\n',...
	    p(i), u(i), v(i),    ddh(i), ddp(i));
  end


  fclose(fid);

end




