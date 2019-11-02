function data=importdata(fname_in)
% for compartibility with matlab

  symlink(fname_in, 'data.link');
  load data.link;
  unlink('data.link');

end
