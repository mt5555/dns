fid=fopen('file000.0000.out')

%time=fscanf(fid,'%e',1);

time=fread(fid,1,'float64');
time

data=fread(fid,4,'float64');
nx=data(1);
ny=data(2);
nz=data(3);
n_var=data(4);

q = fread(fid,nx*ny*nz*n_var,'float64');
q = reshape(q,nx,ny,nz,n_var);
u = squeeze(q(:,:,:,1));
pcolor(u')
fclose(fid)

