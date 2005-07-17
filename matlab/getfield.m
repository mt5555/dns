function [x,y,z,q,time]=getfield(fname);

fidvor=fopen(fname,'r');

  if (fidvor<0) 
    x=0;
    y=0;
    z=0;
    q=zeros([2 2 2]);
    time=-1;
    return; 
  end;
  
  time=fread(fidvor,1,'float64');
  data=fread(fidvor,3,'float64');
  nx=data(1);
  ny=data(2);
  nz=data(3);
  
  x=fread(fidvor,nx,'float64');
  y=fread(fidvor,ny,'float64');
  z=fread(fidvor,nz,'float64');
  
  if (nargout==3) return; end;

  q = fread(fidvor,nx*ny*nz,'float64');
  tmp = fread(fidvor,1,'float64');
  tmp=size(tmp);
  if (tmp(1)~=0) 
    disp('Error reading input file...')
  end
  fclose(fidvor);
  q = reshape(q,nx,ny,nz);

  return
  