%
%########################################################################
%#  plot of DNS spectrum output file
%########################################################################
%
clear
movie=1;       % plot all the spectrum, pausing between each one
movie_plot=0;  % create .ps out of all the spectra
endian = 'l';  % default endian-ness of spec file
CK_orig=1.613/(2*pi)^2;
decay_scale=0;   % special scaling for 2048^3 decaying run
tsave=[];
mu=1;

sc_count=0;
for sc=1:18
sc_count=sc_count+1;

%sc=2;
timename='0000.7019';
%timename='0000.7536';
basename='/scratch1/taylorm/decay2048/decay2048-new.';

namesc=sprintf('%i5',sc+10000); namesc=namesc(2:5);
name=[basename,'-sc',namesc,'_',timename];

% open the 'select' file and look for entry sc:
nameselect=[basename,timename,'.select'];
fid=fopen(nameselect);
for i=1:sc;
  data=fscanf(fid,'%f %f %f %f %f %f %f %f',8);
  il=data(7); jl=data(8);  S=data(6); eps_l=data(5);
  data=fscanf(fid,'%f %f %f',[3,3]); 
  data=data';

  %[data(il,jl),S/2]
end
fclose(fid);
cospec_ind=[0 1 2; -1 0 3; -2 -3 0; ];
cospec_ij = cospec_ind(il,jl);
cospec_scale = S*eps_l^(1/3);
%[il,jl,cospec_ij]

tmp=il; il=jl; jl=tmp;

spec_r_save=[];
spec_r_save_fac3=[];

% note: endianopen() will notwork with .spec files
%because first number is not necessaryly an integer 
fid=fopen([name,'.spec'],'r',endian);
fidco=endianopen([name,'.cospec'],'r');  



time=fread(fid,1,'float64')


CK=CK_orig;
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
  
  ke=sum(spec_r');            % spec_r = .5u_k**2
  ke_diss=mu*2*sum(knum.^2 * (2*pi)^2 .* spec_r')  ; 
  L11 = ((3*pi)/(4*ke)) * sum(spec_r(2:n_r)'./(knum(2:n_r)*2*pi));  

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
      loglog53(n_r,spec_r,stitle,1,4);
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
      %spec_r=spec_r./fac3;
      %loglog53(n_r-1,spec_r,stitle,CK); hold on;
      hold off;
      
      
      % longitudinal spectraum
      figure(4)
      subplot(2,1,1);
      loglog53(n_x,spec_ux,' ',CK*18/55);     hold on;
      loglog53(n_y,spec_vy,' ',CK*18/55);     hold on;
      loglog53(n_z,spec_wz,'longitudinal 1D spectrum',CK*18/55);     hold on;
      hold off;
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
  
  
  
  
  
  


  % now read the cospectrum:
  if (fidco>-1) 
    disp('reading co-spectrum')
    ncox=fread(fidco,1,'float64'); 
    ncoy=fread(fidco,1,'float64'); 
    ncoz=fread(fidco,1,'float64'); 
    ncor=fread(fidco,1,'float64'); 
    time_co=fread(fidco,1,'float64');
    % uv, uw, vw:
    fread(fidco,ncox,'float64');   % 1D spectrum
    fread(fidco,ncox,'float64'); 
    fread(fidco,ncox,'float64'); 
    fread(fidco,ncox,'float64'); 
    fread(fidco,ncox,'float64'); 
    fread(fidco,ncox,'float64'); 
    fread(fidco,ncox,'float64'); 
    fread(fidco,ncox,'float64'); 
    fread(fidco,ncox,'float64'); 
    uv_r=fread(fidco,ncor,'float64');  % 3D spectrum
    uw_r=fread(fidco,ncor,'float64'); 
    vw_r=fread(fidco,ncor,'float64'); 

    
    if (il==1  & jl == 2)              % S12 dominant
      cospec_1=uv_r/cospec_scale;  
      cospec_2=uw_r/cospec_scale;
      cospec_3=vw_r/cospec_scale;
    elseif (il==1 & jl==3)         % S13 dominant
      cospec_1=uw_r/cospec_scale;
      cospec_2=-uv_r/cospec_scale;  
      cospec_3=-vw_r/cospec_scale;
    elseif (il==2 & jl==3)
      cospec_1=vw_r/cospec_scale;
      cospec_2=uv_r/cospec_scale;  
      cospec_3=uw_r/cospec_scale;
    elseif (il==2 & jl==1) 
      cospec_1=uv_r/cospec_scale;  
      cospec_2=-uw_r/cospec_scale;
      cospec_3=vw_r/cospec_scale;
    elseif (il==3 & jl==1)
      cospec_1=uw_r/cospec_scale;
      cospec_2=vw_r/cospec_scale;
      cospec_3=uv_r/cospec_scale;  
    elseif (il==3 & jl==2)
      cospec_1=vw_r/cospec_scale;
      cospec_2=-uw_r/cospec_scale;
      cospec_3=uv_r/cospec_scale;  
    else
      disp('error: invalid cospec_ij value' );
      return;
    end

%    cospec_1=-cospec_1;
    
    if (sc_count==1) 
      mean_1=cospec_1;
      mean_2=cospec_2;
      mean_3=cospec_3;
    else
      mean_1=mean_1+cospec_1;
      mean_2=mean_2+cospec_2;
      mean_3=mean_3+cospec_3;
    end

    
    figure(3); clf;
    subplot(2,1,1);
    loglog(0:ncor-1,cospec_1,'r'); hold on;
    loglog(0:ncor-1,cospec_2,'g');
    loglog(0:ncor-1,cospec_3,'b');
    loglog(0:ncor-1,.01*(0:ncor-1).^-(7/3));
    axis([1 200 1e-9 1e-3]);
    hold off;
    title(sprintf('eps=%.4f  S_{%i,%i}=%.2f   ',eps_l,il,jl,S));
    
    subplot(2,1,2);
    loglog(0:ncor-1,-cospec_1,'r'); hold on;
    loglog(0:ncor-1,-cospec_2,'g');
    loglog(0:ncor-1,-cospec_3,'b');
    loglog(0:ncor-1,.01*(0:ncor-1).^-(7/3));
    axis([1 200 1e-9 1e-3]);
    hold off;

    figure(4); clf;
    

    %   
    mean_scale=sc_count *( (0:ncor-1).^(-7/3))';

    subplot(2,1,1);
    loglog(0:ncor-1,mean_1./mean_scale,'r'); hold on;
    loglog(0:ncor-1,mean_2./mean_scale,'g');
    loglog(0:ncor-1,mean_3./mean_scale,'b');
    loglog(0:ncor-1,.01*(0:ncor-1).^-(7/3));
    axis([1 128 1e-5 1e-1]);
    hold off;
    title(sprintf('-E12 (red),  E13 (green),  E23 (blue)'));
    
    subplot(2,1,2);
    loglog(0:ncor-1,-mean_1./mean_scale,'r'); hold on;
    loglog(0:ncor-1,-mean_2./mean_scale,'g');
    loglog(0:ncor-1,-mean_3./mean_scale,'b');
    axis([1 128 1e-5 1e-1]);
    hold off;
    
  end
  

  
  % make PS files out of plots:
  if (movie_plot==1)  
    if (1) % ( (2*time-floor(2*time))<.01) | (abs(time-2.0)<.01) )
      disp('making ps files ...' )
      figure(1)
      print ('-dpsc',sprintf('%s_%.2f.ps',name,time))
      if (num_spect>0) 
        figure(2)
        print ('-dpsc',sprintf('%s_%.2f_t.ps',name,time))
      end
      disp('pause')
      pause
    else     
      disp('pause')
      pause
    end
  end
  
  
  time=fread(fid,1,'float64');
  time
  
end
fclose(fid);
fclose(fidco);

%'pause'; pause
end
return



if (length(spec_r_save>1) )
  figure(1); clf;
  loglog53(n_r,spec_r_save,'KE spectrum',CK,1);
  print -djpeg -r72 spec.jpg
  print -depsc -r600 spec.ps
  figure(2); clf;
  %loglog53(n_r,spec_r_save_fac3,'KE / bottleneck-factor',CK);
  %print -djpeg -r72 speck3.jpg
  
  k2=0:n_r-1;
  k2=k2.^2;
  spec=diag(k2) * spec_r_save_fac3(1:n_r,:);
  loglog53(n_r,spec,'Enstrophy spectrum',CK,2);
  print -djpeg -r72 enstrophy.jpg
end

