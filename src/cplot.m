fidu=fopen('test-0-0-0-0000.0000.data');

%
%########################################################################
%#  plotting output file
%########################################################################

%ts=input('time=? ');

range=0:50;

for i=range
  ts=i;
  ts = sprintf('%9.5f',10000+ts);
  ts=ts(2:10);
  ts=['test',ts,'.vor'];
  fidvor=fopen(ts,'r');
  time=fread(fidvor,1,'float64')
  data=fread(fidvor,3,'float64');
  nx=data(1);
  ny=data(2);
  nz=data(3);
  
  x=fread(fidvor,nx,'float64');
  y=fread(fidvor,ny,'float64');
  z=fread(fidvor,nz,'float64');
  
  q = fread(fidvor,nx*ny*nz,'float64');
  fclose(fidvor);
  ts=sprintf('time=%f  %ix%ix%i',time,nx,ny,nz);
  
  q = reshape(q,nx,ny,nz);
  
  vor = squeeze(q(:,:,1));
  figure(2)
  pcolor(x,y,vor')
  title(ts);
  shading interp
  axis square
  
end
return




%
%########################################################################
%#  restart file
%########################################################################
time=fread(fidu,1,'float64');


data=fread(fidu,4,'float64');
nx=data(1);
ny=data(2);
nz=data(3);
n_var=data(4);

q = fread(fidu,nx*ny*nz*n_var,'float64');
fclose(fidu);

disp(sprintf('restart dump:\ntime=%f  dims=%i %i %i %i',time,nx,ny,nz,n_var))
q = reshape(q,nx,ny,nz,n_var);

u = squeeze(q(:,:,1,1));
figure(1)
subplot(2,2,1)
pcolor(u')
shading interp

v = squeeze(q(:,:,1,2));
subplot(2,2,2)
pcolor(v')
shading interp


