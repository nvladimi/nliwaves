% Compute propagation for n0 consecutive seeds starting with seed n1
% using input arguments "fbase", "n0" and "n1".

arg_list = argv ();

fbase = arg_list{1};
s1 = arg_list{2};
s2 = arg_list{3};

n1 = str2num(s1);
n0 = str2num(s2);

for seed=n1:n1+n0-1

  giaflucs(fbase, seed);

end

