%
%########################################################################
%#  plot of DNS spectrum output file
%########################################################################
%
CK_orig=1.0;
decay_scale=0;

% name of the file, without the ".spec" extension
name='cj20000.0000';
namedir='/ccs/taylorm/dns/src/';

%name='temp0000.0000';
%namedir='/ccs/taylorm/dns/src/';


name='decay2048';
namedir='/ccs/scratch/taylorm/decay/';
CK_orig=1.613; decay_scale=1;


% plot all the spectrum:
movie=0;


% save spectrum at these times:
tsave=[0 .41 1.0  1.5  2.0  2.5  3.0 3.5 ];
spec_r_save=[];
spec_r_save_fac3=[];

fid=fopen([namedir,name,'.spec'],'r','l');
fidt=endianopen([namedir,name,'.spect'],'r');
fidt=-1;

time=fread(fid,1,'float64');
num_spect=0;
if (fidt>=0) 
  num_spect=fread(fidt,1,'float64');
  time_t=fread(fidt,1,'float64');
end


j=0;
while (time>=.0 & time<=9999.3)
  j=j+1;
  n_r=fread(fid,1,'float64');
  spec_r=fread(fid,n_r,'float64');

  knum=0:(n_r-1);
  eta=0;
  spec_scale=1; 

  if (decay_scale) 
    % convert to 2pi units:
    mu=3.4424e-6*(2*pi)^2;
    eps = mu*2*sum(knum.^2 * (2*pi)^2 .* spec_r')
    eta = (mu^3 ./ eps).^(.25);
    if (j==1) eps_orig=eps; end;
    
    spec_scale=(2*pi)^2*eps^(-2/3);
    
    % make spectrum dimensionless:
    spec_r=spec_r *spec_scale;
    CK=CK_orig;
    
    
    % dimensional scaling:  (undo the above)
    % and use CK from first spectrum
    spec_r=spec_r*eps^(2/3);
    CK=CK_orig*eps_orig^(2/3);
  end

  knum2=knum;
  knum2(1)=.000001;
  fac3=.5 + atan( 10*log10(knum2'*eta)+12.58 ) / pi;
  fac3 = 1 + .522*fac3;


    
  n_x=fread(fid,1,'float64');
  spec_ux=spec_scale*fread(fid,n_x,'float64');
  spec_vx=spec_scale*fread(fid,n_x,'float64');  
  spec_wx=spec_scale*fread(fid,n_x,'float64');  
  
  n_y=fread(fid,1,'float64');
  spec_uy=spec_scale*fread(fid,n_y,'float64');
  spec_vy=spec_scale*fread(fid,n_y,'float64');
  spec_wy=spec_scale*fread(fid,n_y,'float64');
  
  n_z=fread(fid,1,'float64');
  spec_uz=spec_scale*fread(fid,n_z,'float64');
  spec_vz=spec_scale*fread(fid,n_z,'float64');
  spec_wz=spec_scale*fread(fid,n_z,'float64');  

  i=find( abs(time-tsave)<.0001);
  if (length(i)>=1) 
     tsave(i)
     spec_r_save=[spec_r_save, spec_r];
     spec_r_save_fac3=[spec_r_save_fac3, spec_r./fac3];
  end 

  if (movie==1)  
  figure(1);clf;
  if (n_z==1) 
    %subplot(2,1,1);
    subplot(1,1,1);
    stitle=sprintf('Spectrum t=%8.4f',time);
    loglog53(n_r,spec_r,stitle,1);
%    subplot(2,1,2);
%    loglog53(n_x,spec_ux,' ',1);
%    hold on;
%    loglog53(n_y,spec_vy,' ',1);
%    hold off;
  else
    %spherical wave number
    figure(1)
    subplot(1,1,1);
    stitle=sprintf('Spectrum t=%8.4f',time);
    loglog53(n_r-1,spec_r,stitle,CK); hold on;

    spec_r=spec_r./fac3;
    loglog53(n_r-1,spec_r,stitle,CK); hold off;


    
    % longitudinal spectraum
    figure(2)
    subplot(2,1,1);
    loglog53(n_x,spec_ux,' ',CK*18/55);     hold on;
    loglog53(n_y,spec_vy,' ',CK*18/55);     hold on;
    loglog53(n_z,spec_wz,'longitudinal 1D spectrum',CK*18/55);     hold off;
    
    % transverse spectraum

    subplot(2,1,2);
    loglog53(n_x,spec_uy,' ',CK*18/55);     hold on;
    loglog53(n_x,spec_uz,' ',CK*18/55);     hold on;
    loglog53(n_y,spec_vx,' ',CK*18/55);     hold on;
    loglog53(n_y,spec_vz,' ',CK*18/55);     hold on;
    loglog53(n_z,spec_wx,' ',CK*18/55);     hold on;
    loglog53(n_z,spec_wy,'transverse 1D spectrum',CK*18/55);     
    hold off;
  end
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


  if (movie==1)  
  if ( ( (2*time-floor(2*time))<.01) | (abs(time-.45)<.01) )
    disp('making ps files ...' )
    figure(1)
    print ('-dpsc',sprintf('%s_%.2f.ps',name,time))
    if (num_spect>0) 
      figure(2)
      print ('-dpsc',sprintf('%s_%.2f_t.ps',name,time))
    end
    disp('pause')
    %pause
  else     
    disp('pause')
    %pause
  end
  end
  
  time=fread(fid,1,'float64');
  if (fidt>=0) 
    num_spec=fread(fidt,1,'float64');
    time_t=fread(fidt,1,'float64');
  end
end

fclose(fid);
if (fidt>0) fclose(fidt); end;

%return


if (length(spec_r_save>1) )
figure(1); clf;
loglog53(n_r,spec_r_save,'KE spectrum',CK);
print -djpeg -r72 spec.jpg
figure(2); clf;
loglog53(n_r,spec_r_save_fac3,'KE / bottleneck-factor',CK);
print -djpeg -r72 speck3.jpg
end

