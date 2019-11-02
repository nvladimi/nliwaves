
  IP=1
  fname = 'qr3x3.dat'

  IP=2
  fname = 'qr5x5.dat'

  IP=3
  fname = 'qr7x7.dat'

  A=[0,0,0,0];

  for j=-IP:IP 
    for i=-IP:IP 
      A=[A; 1,i,j, i^2+j^2]; 
    endfor 
  endfor

  A = A(2:end,:);

  [q,r]=qr(A);

  r=r(1:4,1:4);

  Qt = q';

  Ri=inv(r);

  myfile = fopen (fname, 'wt');

  m=(2*IP+1)^2;

  for j=1:m
    for i=1:m
      fprintf(myfile, "%20.12e\n", Qt(j,i));
    endfor
  endfor


  for j=1:4
    for i=1:4
      fprintf(myfile, "%20.12e\n", Ri(j,i));
    endfor
  endfor

  fclose(myfile);
