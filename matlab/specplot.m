%
%########################################################################
%#  plot of DNS spectrum output file
%########################################################################
%
%   call cwrite8(fid,time,1)
%   call cwrite8(fid,1+iwave,1)
%   call cwrite8(fid,spectrum,1+iwave)
%   call cwrite8(fid,1+g_nx/2,1)
%   call cwrite8(fid,spec_x,1+g_nx/2)
%   call cwrite8(fid,1+g_ny/2,1)
%   call cwrite8(fid,spec_y,1+g_ny/2)
%   call cwrite8(fid,1+g_nz/2,1)
%   call cwrite8(fid,spec_z,1+g_nz/2)
%

range=0:50;
%range = 1
%range=0:.5:2;

fid=fopen('test64.spec','r');


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
