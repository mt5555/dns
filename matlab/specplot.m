%
%########################################################################
%#  plot of DNS spectrum output file
%########################################################################
%
clear
movie=1;       % plot all the spectrum, pausing between each one
movie_plot=0;  % create .ps out of all the spectra
endian = 'b';  % default endian-ness of spec file
CK_orig=1.0;
decay_scale=0;
tsave=[];



%name='sc1024A';
%namedir='/ccs/scratch/taylorm/dns/sc1024A/';
%CK_orig=1.613;

name='tmix256C';
namedir='/scratch2/taylorm/tmix256C/archive/';
CK_orig=1.613;
movie_plot=0;
endian='l';

%name='decay2048'; namedir='/ccs/scratch/taylorm/decay/';
%CK_orig=1.613; decay_scale=1;
% save spectrum at these times:
%movie=0; tsave=[0 .41 1.0  1.5  2.0  2.5  3.0 3.5 ];



%name = 'sk128_alpha15/sk128_alpha150000.0000';
%namedir = '/home/skurien/dns/src/';

name = 'rot3d_sto0000.0000';
namedir = '/home/taylorm/ccs/dns/src/rot3d/';

%name = 'Rot10000.0000';
%namedir = '/ccs/wingate/Rotation/Rot1/';







spec_r_save=[];
spec_r_save_fac3=[];

% note: endianopen() will notwork with .spec files
%because first number is not necessaryly an integer 
fid=fopen([namedir,name,'.spec'],'r',endian);
fidt=endianopen([namedir,name,'.spect'],'r');
fidp=endianopen([namedir,name,'.pspec'],'r');  
%fidt=-1;
%fidp=-1;

time=fread(fid,1,'float64');
num_spect=0;
if (fidt>=0) 
  num_spect=fread(fidt,1,'float64');
  time_t=fread(fidt,1,'float64');
end


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
    loglog53(n_r-1,spec_r,stitle,CK); hold on;


    
    % longitudinal spectraum
    figure(4)
    subplot(2,1,1);
    loglog53(n_x,spec_ux,' ',CK*18/55);     hold on;
    loglog53(n_y,spec_vy,' ',CK*18/55);     hold on;
    loglog53(n_z,spec_wz,'longitudinal 1D spectrum',CK*18/55);     hold on;
    
    % transverse spectraum

    subplot(2,1,2);
    loglog53(n_x,spec_uy,' ',CK*18/55);     hold on;
    loglog53(n_x,spec_uz,' ',CK*18/55);     hold on;
    loglog53(n_y,spec_vx,' ',CK*18/55);     hold on;
    loglog53(n_y,spec_vz,' ',CK*18/55);     hold on;
    loglog53(n_z,spec_wx,' ',CK*18/55);     hold on;
    loglog53(n_z,spec_wy,'transverse 1D spectrum',CK*18/55);     
   hold on;
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


  % now read the passive scalar spectrum:
  if (fidp>-1) 
     disp('reading passive scalar spectrum')
     npassive=fread(fidp,1,'float64'); 
     time_p=fread(fidp,1,'float64');
     figure(5); clf; subplot(1,1,1)

     np_r=fread(fidp,1,'float64');
     for np=1:npassive  %! 1 ! this should be 1:npassive   after we fix data! 
        pspec_r(:,np)=fread(fidp,np_r,'float64');
     end
     ts=sprintf('passive scalars t=%f',time);
     loglog53(np_r,pspec_r,ts,1.0,3); 

     np_x=fread(fidp,1,'float64');
     for np=1:npassive    
        pspec_x(:,np)=fread(fidp,np_x,'float64');
     end

     np_y=fread(fidp,1,'float64');
     for np=1:npassive    
        pspec_y=fread(fidp,np_y,'float64');
     end
     
     np_z=fread(fidp,1,'float64');
     for np=1:npassive    
        pspec_z=fread(fidp,np_z,'float64');
     end

     np=10;
     x=(0:np_r-1)';
     c2=sum(pspec_r(2:np_r,np));     
     cx2=sum((x.^2).*pspec_r(:,np))/3;     % 3<c,1*c,1>
     cxx2=sum((x.^4).*pspec_r(:,np))/3;    % 3<c,11*c,11> 
     %[c2,cx2,cxx2]
     Gc=c2*(cxx2)/(cx2*cx2);

     disp(sprintf('t=%f  <c> Gc  %f  %f  %f  %f ',time,pspec_r(1,np),Gc))

     hold on;
     x=1:1000;
     xt=1e5*((.001+x).^-5);
     plot(x,xt);
%     xt=1e5*((.001+x).^-7);
%     plot(x,xt);
     hold off; 
% $$$      c2=sum(xt);
% $$$      cx2=sum((x.^2).*xt);
% $$$      cxx2=sum((x.^4).*xt);
% $$$      c2*(cxx2/3)/(cx2*cx2/9)

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
    if (fidp>-1) 
      figure(5)
      print ('-djpeg','-r72',sprintf('pspec%.2f.jpg',time))
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

