function[ux,dir_max] = read_ux(fname)

fid=endianopen(fname,'r')

fname

ux = fread(fid,9,'float64');
ux = reshape(ux,3,3);
ux = ux'
     dir_max = round(fread(fid,9,'float64'))
