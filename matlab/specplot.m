%
%########################################################################
%#  plot of DNS spectrum output file
%########################################################################
%


% name of the file, without the ".spec" extension
name='cj20000.0000';
namedir='/ccs/taylorm/dns/src/';

name='temp0000.0000';
namedir='../src/';

name='decay2048';
namedir='/ccs/scratch/taylorm/decay/';

% save spectrum at these times:
tsave=[0 .41 .80  1.2 1.6 2.0 2.44 2.45 2.6 ];
spec_r_save=[];

fid=endianopen([namedir,name,'.spec'],'r');
fidt=endianopen([namedir,name,'.spect'],'r');
fidt=-1;

time=fread(fid,1,'float64');
num_spect=0;
if (fidt>=0) 
  num_spect=fread(fidt,1,'float64');
  time_t=fread(fidt,1,'float64');
end



j=0;
while (time>=0 & time<=9999.3)
  j=j+1;
  n_r=fread(fid,1,'float64');
  spec_r=fread(fid,n_r,'float64');
  
  if (max(spec_r)<.01) spec_r=spec_r*(2*pi)^2; end;    
    
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

  i=find( abs(time-tsave)<.0001);
  if (length(i)>=1) 
     tsave(i)
     spec_r_save=[spec_r_save, spec_r];
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
  else
    %spherical wave number
    figure(1)
    subplot(1,1,1);
    loglog53(n_r-1,spec_r,time);

    
    % longitudinal spectraum
    figure(2)
    subplot(2,1,1);
    loglog53(n_x,spec_ux,time);     hold on;
    loglog53(n_y,spec_vy,time);     hold on;
    loglog53(n_z,spec_wz,time,'longitudinal 1D spectrum');     hold off;
    
    % transverse spectraum
    subplot(2,1,2);
    loglog53(n_x,spec_uy,time);     hold on;
    loglog53(n_x,spec_uz,time);     hold on;
    loglog53(n_y,spec_vx,time);     hold on;
    loglog53(n_y,spec_vz,time);     hold on;
    loglog53(n_z,spec_wx,time);     hold on;
    loglog53(n_z,spec_wy,time,'transverse 1D spectrum');     
    hold off;
  end


  
  %
  % NOW read the transfer spectrum
  %
  if (num_spect==4) % NS transfer spectrum
    n_r=fread(fidt,1,'float64');
    spec_tot=fread(fidt,n_r,'float64');
    spec_model=0*spec_tot;
    
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

  
  

  if (num_spect>0)
  %
  % transfer spectrum
  %
  figure(3);
  subplot(2,1,1)
  x=0:n_r-1;
  %semilogx(x,spec_transfer,'k',x,spec_diff,'r',x,spec_f,'b');
  semilogx(x,spec_transfer,'k',x,spec_diff,'r',x,spec_model,'c',x,spec_tot,'y');
  title(sprintf('T_k (black)    D_k (red)    M_k(cyan)   dEk/dt(yellow)      time = %f ',time));

  subplot(2,1,2)
  %semilogx(x,spec_transfer+spec_diff+spec_f,x,spec_tot,'o');
  flux=0*spec_transfer;
  for i=1:length(spec_transfer)
     flux(i)=-sum(spec_transfer(1:i)); 
  end      
  semilogx(x,flux); 
  grid;
  title(sprintf('E Flux'));
  end


  
%  if ( ( (2*time-floor(2*time))<.01) | (abs(time-.4020)<.01) )
  if ( time>2.4 )
    disp('making ps files ...' )
    figure(1)
%    print ('-dpsc',sprintf('%s_%.2f.ps',name,time))
    if (num_spect>0) 
      figure(2)
%      print ('-dpsc',sprintf('%s_%.2f_t.ps',name,time))
    end
%    disp('pause')
%    pause
  end     
  

  
  time=fread(fid,1,'float64');
  if (fidt>=0) 
    num_spec=fread(fidt,1,'float64');
    time_t=fread(fidt,1,'float64');
  end
end

fclose(fid);
if (fidt>0) fclose(fidt); end;

if (length(spec_r_save>1) )
figure(1); clf;
loglog53(n_r,spec_r_save,.43,'KE spectrum');
print -djpeg -r72 spec.jpg
end

