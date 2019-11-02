
function guide_nk_shells
%
%  function guide_nk_shells
%
%  Parameters are edited inside.
%
%  Reads *.nk files, computes forcing, source term, and potential
%  for visualizing q-flux.  Data is saved in vtk format for paraview.
%  
%  Subroutine nk_shells prints on the screen source term integrated
%  over pumping and damping shells.  
%  


  fbase='n0800';      N = 256;

 
  dt = 1;


  fk  = guide_forcing(2, N/2);


  fid_rates=fopen([fbase, ".shells_rates.txt"], "wt");
  fid_nk=fopen([fbase, ".shells_nk.txt"], "wt");

  fprintf(fid_rates, "# 1.time  2.n  3.dn_dealias  4.pumping  5.damping  6.inner  7.total  8.k=1\n\n");

  fprintf(fid_nk, "# 1.time  2.n  3.dn_dealias  4.pumping  5.damping  6.inner  7.total  8.k=1\n\n");

  

  %-- loop over files --  

  for k = 0:2000

    num = num2str(k,'%04d');
    %num='avg';

    fname = [fbase, '.psi.', num];

    f=read_psi(fname, N); 
    q=fft2(f)/N/N; 
    nk=q.*conj(q); 

    n0=nk(1,1); 
    nk(1,1)=0;
    n=sum(sum(nk));

    nk = fftshift(nk);
    nk = nk(N/4+1:3*N/4, N/4+1:3*N/4);
    dn = sum(sum(nk)) - n;

    qk = nk.*fk;

    %-- rates --

    nksh = nk_shells(qk, N/2);  
    nksh = [k*dt, n+n0, dn, nksh];
 
    fprintf(fid_rates, "%12.5e  ", nksh);
    fprintf(fid_rates, "\n");

    %-- nk --

    nksh = nk_shells(nk, N/2);  
    nksh = [k*dt, n+n0, dn, nksh];
 
    fprintf(fid_nk, "%12.5e  ", nksh);
    fprintf(fid_nk, "\n");

  end

  fclose(fid_nk);
  fclose(fid_rates);

end



%--------------------------------------------------

function q = nk_shells(f,N)
% N is number of dealiased modes.
% Pumping shell is defined as [kl,kr] including boundaries.

  s = N/128;

  n  = 128*s;
  kl =  28*s;
  kr =  32*s;

 
  %-- grid --

  k = [-n/2:n/2-1]; 

  [kx, ky] = meshgrid (k, k);

  k = sqrt(kx.^2 + ky.^2);

  %-- find total --

  q3 = sum(sum(f));

  %-- erase damping shell --

  ind = find (k>kr);
  f(ind) = 0;

  q2 = sum(sum(f));


  %-- erase pumping shell --

  ind = find (k>=kl);
  f(ind) = 0;

  q1 = sum(sum(f));


  %-- one-cell belt shell --

  ind = find (k<1.5);
  q4 = sum(sum(f(ind)));

  %-- pumping  damping  inner  total square--

  q = [q2-q1, q3-q2, q1, q3, q4]; 

end

%--------------------------------------------------







