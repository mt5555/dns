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
tsave=[];
mu=1;

sc_count=0;
%for sc=1:18
%for sc=1111:1111   
for sc=-1:-1        % regular file, dont use a _sc postfix
sc_count=sc_count+1;

%timename='0000.7019';
%timename='0000.7536';
%basename='/scratch1/taylorm/decay2048/specold/decay2048-new.';
%basename='/scratch1/taylorm/decay2048/decay2048-new.';

%timename='0002.0000';
%basename='/scratch2/taylorm/sc1024A/sc1024A';

%timename='0000.7019';
%basename='/scratch1/taylorm/decay2048/decay2048-new.';

%timename='0000.0000';
%basename='/scratch1/taylorm/shear/shear';

timename='0000.0000';
basename='~/projects/shear/livescu/S=1.275/fort.298';


namesc=sprintf('%i5',sc+10000); namesc=namesc(2:5);
if (sc>=0) 
   name=[basename,'-sc',namesc,'_',timename];
else
   name=[basename,timename];
end

% open the 'select' file and look for entry sc:
nameselect=[basename,timename,'.select'];
il=1; jl=2; S=1; eps_l=1;  % default values if 'select' file doesnt exist
fid=fopen(nameselect);
for i=1:sc;
  [data,n]=fscanf(fid,'%f %f %f %f %f %f %f
 %f',8);
  if (n==8)
    il=data(7); jl=data(8);  S=data(6); eps_l=data(5);
    data=fscanf(fid,'%f %f %f',[3,3]); 
    data=data';
    [il,jl,data(il,jl),S/2]
  end
end
if (fid>=0) fclose(fid); end;


cospec_ind=[0 1 2; -1 0 3; -2 -3 0; ];
cospec_ij = cospec_ind(il,jl);
cospec_scale = S*eps_l^(1/3);
%[il,jl,cospec_ij]

% swap il, jl?
% tmp=il; il=jl; jl=tmp;


% note: endianopen() will notwork with .spec files
%because first number is not necessaryly an integer 
fid=fopen([name,'.spec'],'r',endian);
fidco=endianopen([name,'.cospec'],'r');  

time=1.0;
if (fid>=0) 
  time=fread(fid,1,'float64')
end


CK=CK_orig;
j=0;
while (time>=.0 & time<=9999.3)
  j=j+1;
  
  if (fid>=0) 
    n_r=fread(fid,1,'float64');
    spec_r=fread(fid,n_r,'float64');
    
    knum=0:(n_r-1);
    eta=0;
    spec_scale=1; 
    
    
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
        figure(2)
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
    
  end
  
  % now read the cospectrum:
  if (fidco>-1) 
    disp('reading co-spectrum')
    ncox=fread(fidco,1,'float64'); 
    ncoy=fread(fidco,1,'float64'); 
    ncoz=fread(fidco,1,'float64'); 
    ncor=fread(fidco,1,'float64'); 
    time_co=fread(fidco,1,'float64');

    ncox=ncox+1;
    ncoy=ncoy+1;
    ncoz=ncoz+1;
    % uv, uw, vw:
    uv_x=fread(fidco,ncox,'float64');   % 1D spectrum
    uw_x=fread(fidco,ncox,'float64'); 
    vw_x=fread(fidco,ncox,'float64'); 
    fread(fidco,ncoy,'float64'); 
    fread(fidco,ncoy,'float64'); 
    fread(fidco,ncoy,'float64'); 
    fread(fidco,ncoz,'float64'); 
    fread(fidco,ncoz,'float64'); 
    fread(fidco,ncoz,'float64'); 
    uv_r=fread(fidco,ncor,'float64');  % 3D spectrum
    uw_r=fread(fidco,ncor,'float64'); 
    vw_r=fread(fidco,ncor,'float64'); 

    
    if (il==1  & jl == 2)              % S12 dominant
      cospec_1=uv_r/cospec_scale;  
      cospec_2=uw_r/cospec_scale;
      cospec_3=vw_r/cospec_scale;
    elseif (il==1 & jl==3)         % S13 dominant
      ! u2 -> -u3   u3 -> u2   u1,3 -> u1,2
      ! E12 -> -E13    E13 -> E12,  E23 -> -E23
      cospec_1=uw_r/cospec_scale;
      cospec_2=-uv_r/cospec_scale;  
      cospec_3=-vw_r/cospec_scale;
    elseif (il==2 & jl==3)
      ! u1 -> u2   u2 -> u3   u3 -> u1  
      cospec_1=vw_r/cospec_scale;
      cospec_2=uv_r/cospec_scale;  
      cospec_3=uw_r/cospec_scale;
    elseif (il==2 & jl==1) 
      ! u2 -> -u1   u1 -> u2   u2,1 -> -u1,2
      ! E12 -> -E12    E13 -> E23,  E23 -> -E13   
      ! add extra sign change because sign of shear changed:
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
    axis([1 200 1e-9 1e-0]);
    %axis([1 200 1e-15 1e-3]);
    hold off;
    title(sprintf('eps=%.4f  S_{%i,%i}=%.2f   ',eps_l,il,jl,S));
    
    subplot(2,1,2);
    loglog(0:ncor-1,-cospec_1,'r'); hold on;
    loglog(0:ncor-1,-cospec_2,'g');
    loglog(0:ncor-1,-cospec_3,'b');
    loglog(0:ncor-1,.01*(0:ncor-1).^-(7/3));
    axis([1 200 1e-9 1e-0]);
    %axis([1 200 1e-15 1e-3]);
    hold off;

    figure(3)
    orient tall;
    pname=sprintf('cospec%i.ps',sc_count);
    print ('-dpsc',pname);
    
    
    figure(4); clf;
    %   
    mean_scale=sc_count;%  *( (0:ncor-1).^(-7/3))';

    subplot(2,1,1);
    loglog(0:ncor-1,mean_1./mean_scale,'r'); hold on;
    loglog(0:ncor-1,mean_2./mean_scale,'g');
    loglog(0:ncor-1,mean_3./mean_scale,'b');
    loglog(0:ncor-1,.01*(0:ncor-1).^-(7/3),'k');
%    axis([1 128 1e-9 1e-3]);
    hold off;
    title(sprintf('E12 (red),  E13 (green),  E23 (blue)'));
    
    subplot(2,1,2);
    loglog(0:ncor-1,-mean_1./mean_scale,'r'); hold on;
    loglog(0:ncor-1,-mean_2./mean_scale,'g');
    loglog(0:ncor-1,-mean_3./mean_scale,'b');
     loglog(0:ncor-1,.01*(0:ncor-1).^-(7/3),'k');
%    axis([1 128 1e-9 1e-3]);
    hold off;
    
    figure(5);clf;
    loglog(0:ncor-1,abs(mean_1./mean_scale),'r'); hold on;
    loglog(0:ncor-1,abs(mean_2./mean_scale),'g');
    loglog(0:ncor-1,abs(mean_3./mean_scale),'b');
    hold off;
    title(sprintf('|E12| (red),  |E13| (green),  |E23| (blue)'));
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
  
  if (fid>=0)   
    time=fread(fid,1,'float64');
  else
    time=-1;
  end
  time
  
end
if(fid>=0) fclose(fid); end;
fclose(fidco);
end



figure(4)
orient tall
print -dpsc cospecave.ps


