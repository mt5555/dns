fidu=fopen('data000.5000.out');
fidvor=fopen('vor001.9000.out');


%time=fscanf(fidu,'%e',1);

time=fread(fidu,1,'float64');
time

data=fread(fidu,4,'float64');
nx=data(1);
ny=data(2);
nz=data(3);
n_var=data(4);

q = fread(fidu,nx*ny*nz*n_var,'float64');
fclose(fidu);


q = reshape(q,nx,ny,nz,n_var);

u = squeeze(q(:,:,:,1));
figure(1)
subplot(2,2,1)
pcolor(u')
shading interp

v = squeeze(q(:,:,:,2));
subplot(2,2,2)
pcolor(v')
shading interp



time=fread(fidvor,1,'float64');
time

data=fread(fidvor,4,'float64');
nx=data(1);
ny=data(2);
nz=data(3);
n_var=data(4);

q = fread(fidvor,nx*ny*nz*n_var,'float64');
fclose(fidvor);

q = reshape(q,nx,ny,nz,n_var);
vor = squeeze(q(:,:,:,1));
figure(2)
pcolor(vor')
shading interp



