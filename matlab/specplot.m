%
%########################################################################
%#  plot of DNS spectrum output file
%########################################################################
%



%fid=fopen('test32.spec','r','l');
%fid=fopen('iso12_256_200.spec','r','b');

%fid=fopen('../src/impulse/kh4.spec','r','l');
%fid=fopen('../src/kh/khK.spec','r','l');
%fid=fopen('../src/sht/rung_5000_0000.0000.spec','r','l');

fid=fopen('../src/temp0000.0000.spec');
fidt=fopen('../src/temp0000.0000.spect');



time=fread(fid,1,'float64');
num_spec=fread(fidt,1,'float64');
time_t=fread(fidt,1,'float64');


j=0;
while (time>=0 & time<=9999.3)
  j=j+1;
  n_r=fread(fid,1,'float64');
  spec_r=fread(fid,n_r,'float64');
  n_x=fread(fid,1,'float64');
  spec_ux=fread(fid,n_x,'float64');  
  spec_vx=fread(fid,n_x,'float64');  
  spec_wx=fread(fid,n_x,'float64');  
  
  n_y=fread(fid,1,'float64');
  spec_uy=fread(fid,n_y,'float64');
  spec_vy=fread(fid,n_y,'float64');
  spec_wy=fread(fid,n_y,'float64');
  
  n_z=fread(fid,1,'float64');
  spec_uz=fread(fid,n_z,'float64');
  spec_vz=fread(fid,n_z,'float64');
  spec_wz=fread(fid,n_z,'float64');  
  
  %
  % NOW read the transfer spectrum
  %
  n_r=fread(fidt,1,'float64');
  spec_tot=fread(fidt,n_r,'float64');
   
  time_terms=fread(fidt,1,'float64');  
  n_r=fread(fidt,1,'float64');
  spec_transfer=fread(fidt,n_r,'float64');

  time_terms=fread(fidt,1,'float64');  
  n_r=fread(fidt,1,'float64');
  spec_diff=fread(fidt,n_r,'float64');

  time_terms=fread(fidt,1,'float64');  
  n_r=fread(fidt,1,'float64');
  spec_f=fread(fidt,n_r,'float64');

  % if num_spec>4, read the rest of the spectrums:
  if (num_spec>4) 
    for i=1:num_spec-4
      time_terms=fread(fidt,1,'float64');  
      n_r=fread(fidt,1,'float64');
      spec_dummy=fread(fidt,n_r,'float64');
    end
  end


  
  
  
  if (n_z==1) 
    figure(4);
    subplot(2,1,1);
    loglog53(n_r,spec_r,time);
    subplot(2,1,2);
    loglog53(n_x,spec_ux,time);
    hold on;
    loglog53(n_y,spec_vy,time);
    hold off;
    %print -depsc spec.ps    
    %pause
  else
    figure(4);
    clf
    %spherical wave number
    subplot(3,1,1);
    loglog53(n_r-1,spec_r,time);
    
    % longitudinal spectraum
    subplot(3,1,2);
    loglog53(n_x,spec_ux,time);     hold on;
    loglog53(n_y,spec_vy,time);     hold on;
    loglog53(n_z,spec_wz,time,'longitudinal 1D spectrum');     hold off;
    
    % transverse spectraum
    subplot(3,1,3);
    loglog53(n_x,spec_uy,time);     hold on;
    loglog53(n_x,spec_uz,time);     hold on;
    loglog53(n_y,spec_vx,time);     hold on;
    loglog53(n_y,spec_vz,time);     hold on;
    loglog53(n_z,spec_wx,time);     hold on;
    loglog53(n_z,spec_wy,time,'transverse 1D spectrum');     
    hold off;
    
    
  end

  figure(1);
  subplot(2,1,1)
  x=0:n_r-1;
  semilogx(x,spec_transfer,'k',x,spec_diff,'r',x,spec_f,'b');
  title(sprintf('time = %f ',time));
  subplot(2,1,2)
  semilogx(x,spec_transfer+spec_diff+spec_f,x,spec_tot,'o');
  
  'pause...'
  pause

  
  time=fread(fid,1,'float64');
  num_spec=fread(fidt,1,'float64');
  time_t=fread(fidt,1,'float64');
end

fclose(fid);
fclose(fidt);
