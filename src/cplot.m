fidu=fopen('test-0-0-0-0000.0000.data');
%fidvor=fopen('test0000.5000.vor','r');
fidvor=fopen('test0020.0000.vor','r');

%
%########################################################################
%#  plotting output file
%########################################################################
time=fread(fidvor,1,'float64')

data=fread(fidvor,3,'float64');
nx=data(1);
ny=data(2);
nz=data(3);

x=fread(fidvor,nx,'float64');
y=fread(fidvor,ny,'float64');
z=fread(fidvor,nz,'float64');

q = fread(fidvor,nx*ny*nz,'float64');
temp=fread(fidvor,1,'float64')
fclose(fidvor);
disp(sprintf('output file:\ntime=%f  dims=%i %i %i',time,nx,ny,nz))

q = reshape(q,nx,ny,nz);

vor = squeeze(q(:,:,1));
figure(2)
pcolor(x,y,vor')
shading interp
axis square






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


