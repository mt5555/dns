function [ndelta,ndir,r_val,w2w2,s2s2,w2s2] = readisow2s2(fname)


fid=endianopen(fname,'r');

time=fread(fid,1,'float64');
ndelta=fread(fid,1,'float64');
ndir  =fread(fid,1,'float64');
ncor  =fread(fid,1,'float64');
r0means = fread(fid,ncor,'float64');

ncor
r0means

r_val=fread(fid,[ndelta,ndir],'float64');
if (ncor==3) 
   w2w2=fread(fid,[ndelta,ndir],'float64');
   s2s2=fread(fid,[ndelta,ndir],'float64');
   w2s2=fread(fid,[ndelta,ndir],'float64');
end



r_val = [ zeros([1,ndir]) ; r_val];

r0 = r0means(1)*ones([1,ndir]);
w2w2 = [ r0  ; w2w2];

r0 = r0means(2)*ones([1,ndir]);
s2s2 = [ r0  ; s2s2];

r0 = r0means(3)*ones([1,ndir]);
w2s2 = [ r0  ; w2s2];

return
