%
%########################################################################
%#  plot of DNS spectrum output file
%########################################################################
%



%fid=fopen('test32.spec','r','l');
%fid=fopen('iso12_256_200.spec','r','b');

%fid=fopen('../src/impulse/kh4.spec','r','l');
%fid=fopen('../src/kh/khK.spec','r','l');
fid=fopen('../src/sht/rung_5000_0000.0000.spec','r','l');


time=fread(fid,1,'float64');
j=0;
while (time>=0 & time<=9999.3)
  j=j+1;
  n_r=fread(fid,1,'float64');
  spec_r=fread(fid,n_r,'float64');
  n_x=fread(fid,1,'float64');
  spec_x=fread(fid,n_x,'float64');
  n_y=fread(fid,1,'float64');
  spec_y=fread(fid,n_y,'float64');
  n_z=fread(fid,1,'float64');
  spec_z=fread(fid,n_z,'float64');

  if (n_z==1) 
    figure(4);
    subplot(3,1,1);
    loglog53(n_r,spec_r,time);
    subplot(3,1,2);
    loglog53(n_x,spec_x,time);
    subplot(3,1,3);
    loglog53(n_y,spec_y,time);
    %print -depsc spec.ps    
    %pause
  else
    figure(4);
    subplot(2,2,1);
    loglog53(n_r-1,spec_r,time);
    subplot(2,2,2);
    loglog53(n_x,spec_x,time);
    subplot(2,2,3);
    loglog53(n_y,spec_y,time);
    subplot(2,2,4);
    loglog53(n_z,spec_z,time);
    %M(j)=getframe(gcf);
  end
  time=fread(fid,1,'float64');
end

fclose(fid);
