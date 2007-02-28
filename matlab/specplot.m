%
%########################################################################
%#  plot of DNS spectrum output file
%########################################################################
%
clear
movie=1;       % plot all the spectrum, pausing between each one
movie_plot=0;  % create .ps out of all the spectra
movie_times=[];% create .ps files out of only these times
               % if [], create .ps file out of every plot
endian = 'l';  % default endian-ness of spec file
CK_orig=1.0;
decay_scale=0;   % special scaling for 2048^3 decaying run
tsave=[];
mu = 1;

%name='decay2048-ave64.-new.0000.4026';
%name='decay2048-ave32.-new.0000.4026';
%namedir='/home/taylorm/rider/';
%movie_plot=0; CK_orig=1.613; decay_scale=1; endian='l';


%  name='decay2048'; namedir='/home/mataylo/codes/dns_data/decay/';
%  CK_orig=1.613; decay_scale=1; endian='l';
%  % save spectrum at these times:
%  movie_plot=0; movie=1; 
%  tsave=[0 .41 1.0  1.5  2.0  2.5  3.0 3.5 ];
%  tsave=[0 .41 .7 1.0  1.3  1.6 ];


%name = 'qg64hyper_all';
%namedir = '/nh/u/skurien/projects/pv/data_analysis/lowforc/low4/qg64/sto_high_4/hyper_nu/bous500/';


%name = 'qg256hyper_all';
%namedir = '/nh/u/skurien/projects/pv/data_analysis/lowforc/low4/qg256/bous1000/hyper_nu2.5/';


%name = 'qg256hyper0000.0000';
%namedir = '/nh/u/skurien/projects/pv/data_analysis/lowforc/low4/qg256/bous1000/hyper_nu15/';

%name = 'qg256hyper_all';
%namedir = '/nh/u/skurien/projects/pv/data_analysis/lowforc/low4/qg256/bous500/hyper_nu2.5/';

name = 'qg256hyper_all';
namedir = '/nh/u/skurien/projects/pv/data_analysis/lowforc/low4/qg256/bous2000/';

%name = 'sc1024A';
%namedir = '/nh/xraid/skurien/dns_data/sc1024a/';

%name = 'n256_f2000n5_all';
%namedir = '/nh/u/skurien/projects/pv/data_analysis/lowforc/low4/Ro0Fr1/n256/';

%name = 'n256_f5n2000_all';
%namedir = '/nh/u/skurien/projects/pv/data_analysis/lowforc/low4/Ro1Fr0/n256/';



spec_r_save=[];
spec_r_save_fac3=[];

% note: endianopen() will notwork with .spec and .hspec files
%because first number is not necessaryly an integer 
[namedir,name]
fid=fopen([namedir,name,'.spec'],'r',endian);
fidt=endianopen([namedir,name,'.spect'],'r');
fidp=endianopen([namedir,name,'.pspec'],'r');  
fidco=endianopen([namedir,name,'.cospec'],'r');  


fidt=-1;
%fidp=-1;

time=fread(fid,1,'float64')
num_spect=0;
if (fidt>=0) 
  num_spect=fread(fidt,1,'float64');
  time_t=fread(fidt,1,'float64');
end


CK=CK_orig;
j=0;
while (time>=.0 & time<=9000.6)
  j=j+1;
  n_r=fread(fid,1,'float64');
  spec_r=fread(fid,n_r,'float64');
  n_r

  knum=0:(n_r-1);
  eta=0;
  spec_scale=1; 

  if (decay_scale) 
    % convert to 2pi units:
    mu=3.4424e-6*(2*pi)^2;
    eps = mu*2*sum(knum.^2 * (2*pi)^2 .* spec_r');
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
  %n_x
  
  n_y=fread(fid,1,'float64');
  spec_uy=spec_scale*fread(fid,n_y,'float64');
  spec_vy=spec_scale*fread(fid,n_y,'float64');
  spec_wy=spec_scale*fread(fid,n_y,'float64');
  %n_y
  
  n_z=fread(fid,1,'float64');
  spec_uz=spec_scale*fread(fid,n_z,'float64');
  spec_vz=spec_scale*fread(fid,n_z,'float64');
  spec_wz=spec_scale*fread(fid,n_z,'float64');  
  %n_z
  
  i=find( abs(time-tsave)<.0001);
  if (length(i)>=1) 
     tsave(i)
     spec_r_save=[spec_r_save, spec_r];
     spec_r_save_fac3=[spec_r_save_fac3, spec_r./fac3];
  end 

  if (movie==1)  
  figure(1);
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
    stitle=sprintf('Kinetic energy shell-averaged spectrum t=%8.4f',time);
    
    loglog53(n_r-1,spec_r,stitle,CK);hold off;
%   loglog((0:n_r-1),spec_r.*[0:n_r-1]'.^(5/3),'b', 'linewidth',3); hold off;      
    
    %spec_r=spec_r./fac3;
    %loglog53(n_r-1,spec_r,stitle,CK); hold on;
    

    
    % longitudinal spectraum
%    figure(4);title('longitudinal 1D spectrum');
%    subplot(2,1,1);
%    loglog53(n_x,spec_ux,' ',CK*18/55);     hold off;
%    loglog53(n_y,spec_vy,' ',CK*18/55);     hold on;
%    loglog53(n_z,spec_wz,'longitudinal 1D spectrum',CK*18/55);     hold on;
%    hold off;
    
%    % transverse spectraum
%    subplot(2,1,2);
%    loglog53(n_x,spec_uy,' ',CK*18/55);     hold on;
%   loglog53(n_x,spec_uz,' ',CK*18/55);     hold off;
%    loglog53(n_y,spec_vx,' ',CK*18/55);     hold on;
%    loglog53(n_y,spec_vz,' ',CK*18/55);     hold off;
%    loglog53(n_z,spec_wx,' ',CK*18/55);     hold on;
%    loglog53(n_z,spec_wy,'transverse 1D spectrum',CK*18/55);     
%    hold off;
%     title('longitudinal 1D spectra of kinetic energy');

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
     np_r
     pspec_r=[];
     for np=1:npassive 
        pspec_r(:,np)=fread(fidp,np_r,'float64');
        c2(np)=sum(pspec_r(:,np)); 
        % c2_diss = d/dt of .5<c^2>
        c2_diss(np) = mu*2*sum(knum.^2 * (2*pi)^2 .* pspec_r(:,np)')  ; 
        L11c(np) = ((3*pi)/(4*c2(np))) * sum(pspec_r(2:np_r,np)'./ (knum(2:np_r)*2*pi))  ;
     end
     ts=sprintf('shell averaged passive scalar spectrum t=%f',time);
     loglog53(np_r,pspec_r,ts,1.0,3); hold off;
     axis([1 100 1e-6 1]);

     figure(9);hold off
     loglog53(np_r,pspec_r,ts,1.0,3); hold on;
     np_x=fread(fidp,1,'float64');
     for np=1:npassive    
        pspec_x(:,np)=fread(fidp,np_x,'float64');
     end
     ts=sprintf('plane averaged  passive scalars t=%f',time);
     loglog53(np_x, pspec_x, ts, 1.0,3);hold on;    
     
     np_y=fread(fidp,1,'float64');
     for np=1:npassive    
        pspec_y=fread(fidp,np_y,'float64');
     end
     ts=sprintf('plane averaged passive scalars t=%f',time);
     loglog53(np_y, pspec_y, ts, 1.0,3);hold on;
     
     np_z=fread(fidp,1,'float64');
     for np=1:npassive    
        pspec_z=fread(fidp,np_z,'float64');
     end
     ts=sprintf('plane averaged passive scalars t=%f',time);
     loglog53(np_z, pspec_z, ts, 1.0,3);hold off;
     axis([1 100  1e-6 1]);

     np=1;
     % 10  .508
     % 
     % 1   .499905
     x=(0:np_r-1)';
     c2=sum(pspec_r(2:np_r,np));     
     cx2=sum((x.^2).*pspec_r(:,np))/3;     % 3<c,1*c,1>
     cxx2=sum((x.^4).*pspec_r(:,np))/3;    % 3<c,11*c,11> 
     %[c2,cx2,cxx2]
     Gc=c2*(cxx2)/(cx2*cx2);

     disp(sprintf('t=%f  <c> Gc  %f  %f  %f  %f ',time,sqrt(pspec_r(1,np)),Gc))


       % plot the total energy spectrum
  figure(6);
  ts=sprintf('total energy t=%f',time); 
  loglog53(np_r,spec_r+pspec_r,ts,1.0,3);hold off;
  end


  
  
  % now read the cospectrum:
  if (fidco>-1) 
     disp('reading co-spectrum')
     [ncox,count]=fread(fidco,1,'float64'); 
     ncoy=fread(fidco,1,'float64'); 
     ncoz=fread(fidco,1,'float64'); 
     ncor=fread(fidco,1,'float64'); 
     time_co=fread(fidco,1,'float64');
     ncox=ncox+1;
     ncoy=ncoy+1;
     ncoz=ncoz+1;
     % uv, uw, vw:
     fread(fidco,ncox,'float64');   % 1D spectrum
     fread(fidco,ncox,'float64'); 
     fread(fidco,ncox,'float64'); 
     fread(fidco,ncoy,'float64'); 
     fread(fidco,ncoy,'float64'); 
     fread(fidco,ncoy,'float64'); 
     fread(fidco,ncoz,'float64'); 
     fread(fidco,ncoz,'float64'); 
     fread(fidco,ncoz,'float64'); 
     uv_r=fread(fidco,ncor,'float64');  % 3D spectrum
     uw_r=fread(fidco,ncor,'float64'); 
     vw_r=fread(fidco,ncor,'float64'); 

     figure(3); clf;
     subplot(3,1,1);
     loglog(0:ncor-1,uv_r,'r'); hold on;
     loglog(0:ncor-1,uw_r,'g');
     loglog(0:ncor-1,vw_r,'b');
     loglog(0:ncor-1,(0:ncor-1).^-(7/3));
%     axis([1 200 1e-8 1e-2]);
     hold off;

     subplot(3,1,2);
     loglog(0:ncor-1,-uv_r,'r'); hold on;
     loglog(0:ncor-1,-uw_r,'g');
     loglog(0:ncor-1,-vw_r,'b');
     loglog(0:ncor-1,(0:ncor-1).^-(7/3));
%     axis([1 200 1e-8 1e-2]);
     hold off;

     subplot(3,1,3);
     loglog(0:ncor-1,uv_r,'r'); hold on;
     loglog(0:ncor-1,uw_r,'g');
     loglog(0:ncor-1,vw_r,'b');
     loglog(0:ncor-1,-uv_r,'r'); hold on;
     loglog(0:ncor-1,-uw_r,'g');
     loglog(0:ncor-1,-vw_r,'b');
     loglog(0:ncor-1,(0:ncor-1).^-(7/3));
%     axis([1 200 1e-8 1e-2]);
     hold off;
  end



  % make PS files out of plots:
  if (movie_plot==1)  
  if (isempty(movie_times) |  1==max(abs(time-movie_times)<.001) )    
    disp('making ps files ...' )
    figure(1)
    print ('-depsc',sprintf('%s_%.2f.eps',name,time))
    if (num_spect>0) 
      figure(2)
      print ('-depsc',sprintf('%s_%.2f_t.eps',name,time))
    end
    if (fidp>-1) 
      figure(5)
      print ('-depsc',sprintf('pspec%.2f.eps',time))
    end       
    disp('pause')
    pause
  end
  end
  disp('pause')
  pause

  if (fidp>-1)   
    [L11c(1),L11,ke,ke_diss,c2];
    ak(j)=ke_diss*L11/(2*ke^(1.5)/3);
    ac_temp=(c2_diss./c2) .* ((L11c.^2 / ke_diss).^(1/3));
    ac(j)=ac_temp(1);
    time_a(j)=time;
  end
  
  
  time=fread(fid,1,'float64');
  %time
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
loglog53(n_r,spec_r_save,'KE spectrum',CK,5);
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

