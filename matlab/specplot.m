%
%########################################################################
%#  plot of DNS spectrum output file
%########################################################################
%

range=0:50;
%range = 1
%range=0:.5:2;

fid=fopen('test32.spec','r','l');
%fid=fopen('n128.spec','r','b');


time=fread(fid,1,'float64');
j=0;
while (time>=0 & time<=5000)
  j=j+1;
  n_r=fread(fid,1,'float64');
  spec_r=fread(fid,n_r,'float64');
  n_x=fread(fid,1,'float64');
  spec_x=fread(fid,n_x,'float64');
  n_y=fread(fid,1,'float64');
  spec_y=fread(fid,n_y,'float64');
  n_z=fread(fid,1,'float64');
  spec_z=fread(fid,n_z,'float64');

%  figure(1);
%  loglog53(n_r,spec_r,time)
  
  
  figure(4);
  subplot(2,2,1);
  loglog53(n_r,spec_r,time);
  subplot(2,2,2);
  loglog53(n_x,spec_x,time);
  subplot(2,2,3);
  loglog53(n_y,spec_y,time);
  subplot(2,2,4);
  loglog53(n_z,spec_z,time);
  %M(j)=getframe(gcf);
  
  time=fread(fid,1,'float64');
end

fclose(fid);
