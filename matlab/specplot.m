%
%########################################################################
%#  plot of DNS spectrum output file
%########################################################################
%




%fid=fopen('/ccs/scratch/taylorm/dns/iso12_512.spec','r','b');
%fidt=fopen('/ccs/scratch/taylorm/dns/iso12_512.spect','r','b');

fid=fopen('../src/temp0000.0000.spec');
fidt=fopen('../src/temp0000.0000.spect');




time=fread(fid,1,'float64');
num_spect=fread(fidt,1,'float64');
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
  if (num_spect==4) % NS transfer spectrum
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
  end

  if (num_spect==2) % Shallow water transfer spectrum
    n_r=fread(fidt,1,'float64');
    spec_tot=fread(fidt,n_r,'float64');
    spec_f=0*spec_tot;    
    spec_model=0*spec_tot;

    time_terms=fread(fidt,1,'float64');  
    n_r=fread(fidt,1,'float64');
    spec_diff=fread(fidt,n_r,'float64');

    spec_transfer=spec_tot-spec_diff;
  end
  if (num_spect==3) % Shallow water & modeling transfer spectrum
    n_r=fread(fidt,1,'float64');
    spec_tot=fread(fidt,n_r,'float64');
    spec_f=0*spec_tot;    

    time_terms=fread(fidt,1,'float64');  
    n_r=fread(fidt,1,'float64');
    spec_diff=fread(fidt,n_r,'float64');

    time_terms=fread(fidt,1,'float64');  
    n_r=fread(fidt,1,'float64');
    spec_model=fread(fidt,n_r,'float64');

    spec_transfer=spec_tot-spec_diff-spec_model;
  end

  
  
  
  figure(1);clf;
  if (n_z==1) 
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
    orient tall
    print -depsc spec.ps    
    pause

  end

  %
  % transfer spectrum
  %
  figure(2);
  subplot(2,1,1)
  x=0:n_r-1;
  %semilogx(x,spec_transfer,'k',x,spec_diff,'r',x,spec_f,'b');
  semilogx(x,spec_transfer,'k',x,spec_diff,'r',x,spec_model,'y');
  title(sprintf('T_k (black)      D_k (red)        time = %f ',time));
  subplot(2,1,2)
  %semilogx(x,spec_transfer+spec_diff+spec_f,x,spec_tot,'o');
  flux=0*spec_transfer;
  for i=1:length(spec_transfer)
     flux(i)=-sum(spec_transfer(1:i)); 
  end      
  semilogx(x,flux); 
  grid;
  title(sprintf('E Flux'));
  
  'pause...'
  pause

  
  time=fread(fid,1,'float64');
  num_spec=fread(fidt,1,'float64');
  time_t=fread(fidt,1,'float64');
end

fclose(fid);
fclose(fidt);
