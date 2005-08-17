function[ux,dir_max] = read_ux(fname)

fid=fopen(fname,'r')

fname

time = fread(fid,1,'float64')
ux = fread(fid,9,'float64');
ux = reshape(ux,3,3);
format short g
ux = ux';
     dir_max = round(fread(fid,1,'float64'))
