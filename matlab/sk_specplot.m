%
%########################################################################
%#  plot of DNS spectrum output file
%########################################################################
%
CK_orig=1.0;
decay_scale=0;
tsave=[2.5];

% name of the file, without the ".spec" extension
%num = sprintf('%d',t+100)
%name = ['sk128_alpha',num(2:3), '/sk128_alpha',num(2:3),'0000.0000']
%namedir = '/home/skurien/dns/src/';

%name = 'sc1024A0001.9912';
%namedir = '/home2/skurien/helicity_data/sc1024_data/';

%name = 'sk128_hq_0000.0000';
%namedir = '/home2/skurien/helicity_data/helical_forced/';

%name = 'sk128_alpha000000.0000'
%namedir = '/home2/skurien/dns/src/sk128_alpha00/v5e-4/';

%name = 'hel128_hpi4_0000.0000';
%namedir = '/home2/skurien/dns/src/hel128_hpi4/';

%name = 'hel128_hpi2_0000.0000';
%namedir = '/home2/skurien/dns/src/';

name = 'hel256_hpi2_all';
namedir = '/home2/skurien/helicity_data/helical_forced/hel256_hpi2/';
%name = 'skhel512a0000.0000';
%namedir = '/home/mt/data/skhel/';

mu = 2e-4;

%name = 'Rot10000.0000';
%namedir = '/home2/skurien/rotation/Rot1/';

%name = 'sk128_v5e-4_alpha400000.0000';
%namedir = '/home2/skurien/dns/src/sk128_alpha40/';

%name = 'hel480_hpi2_0000.1000';
%namedir = '/home2/skurien/helicity_data/helical_forced/hel480_hpi2/';
%mu = 2e-4;

%name = 'hel256_h0_0000.8000';
%namedir = '/home2/skurien/helicity_data/helical_forced/hel256_h0/';
%mu = 2e-4;

%namedir = '/home2/skurien/dns/src/';
%name = '2Dhypo1e3_1024_0000.0000';
%mu = 2e-4;

% plot all the spectrum:
movie=1; % 1 to plot all the spectra

spec_r_save=[];
spec_r_save_fac3=[];
pname = [strrep(name,'_','-'),'.spec'];

disp([namedir,name,'.spec']);
fid=fopen([namedir,name,'.spec'],'r','l');
fidt=endianopen([namedir,name,'.spect'],'r');
%fidt=-1;

time=fread(fid,1,'float64');
time
num_spect=0;
if (fidt>=0) 
  num_spect=fread(fidt,1,'float64');
  time_t=fread(fidt,1,'float64');
end
  spec_r = [];
  spec_ave = spec_r;

CK=CK_orig;
j=0;
while (time>=.0 & time<=9999.3)
  n_r=fread(fid,1,'float64');
  spec_r=fread(fid,n_r,'float64');
%if (j==0)
if (time < 1.5)
     spec_ave =  0*spec_r;  %initialize spec_ave
end
if (time >= 1.5)        %only average times greater than 1 eddy turnover
     j=j+1;                 %count j's for averaging only after t==1.
     spec_ave = spec_ave + spec_r;
end

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
     tsave(i);
     spec_r_save=[spec_r_save, spec_r];
     spec_r_save_fac3=[spec_r_save_fac3, spec_r./fac3];
  end 

  if (movie==1)  
  figure(1);
  if (n_z==1) 
    stitle=sprintf('Spectrum t=%8.4f',time);
    loglog53(n_r,spec_r,stitle,1);
  else
    %sphevrical wave number
    
    figure(1)
   
     stitle=sprintf('Spectrum t=%8.4f',time);
    
    loglog53(n_r-1,spec_r,stitle,CK); hold on;
    spec_r=spec_r./fac3;
    loglog53(n_r-1,spec_r,stitle,CK); hold on;

    
    % longitudinal spectraum
    figure(2)
    loglog53(n_x,spec_ux,' ',CK*18/55);     hold on;
    loglog53(n_y,spec_vy,' ',CK*18/55);     hold on;
    loglog53(n_z,spec_wz,'longitudinal 1D spectrum',CK*18/55);
    hold on;
   
    % transverse spectraum
    figure(3)
    loglog53(n_x,spec_uy,' ',CK*18/55);     hold on;
    loglog53(n_x,spec_uz,' ',CK*18/55);     hold on;
    loglog53(n_y,spec_vx,' ',CK*18/55);     hold on;
    loglog53(n_y,spec_vz,' ',CK*18/55);     hold on;
    loglog53(n_z,spec_wx,' ',CK*18/55);     hold on;
    loglog53(n_z,spec_wy,'transverse 1D spectrum',CK*18/55);     
   hold on;

    %compensated spectrum
   figure(4)
     loglog((0:n_r-1),spec_r.*[0:n_r-1]'.^(5/3),'b'); hold on;
     title('5/3 Compensated spectra');

   
    
  end
  end
   time=fread(fid,1,'float64');
   time
%compute total Energy for snapshot
E = sum(spec_r);
disp(sprintf('Total Energy  = %d',E));

k = [0:n_r-1];
%compute dissipation rate of energy
epsilon = mu*sum(spec_r .* (k'*2*pi).^2);
disp(sprintf('energy dissipation rate = %d',epsilon));

end
fclose(fid);  

%end

% mean energy spectrum
spec_ave = spec_ave/j;

figure(5)
loglog(k,spec_ave,'k--'); hold on;
title('Mean energy spectrum');
ylabel(pname);
xlabel('k')

% compensated mean spectra

figure(6)
loglog(k,spec_ave.*k'.^(5/3),'b'); hold on;
loglog(k,spec_ave.*k'.^(4/3),'r');hold on;
legend('5/3 comp','4/3 comp');
     title('Compensated MEAN energy spectrum');
ylabel(pname);
xlabel('k');



%compute mean energy
E = sum(spec_ave);
disp(sprintf('Mean energy = %d',E));

%compute mean dissipation rate
epsilon = mu*sum(spec_ave .* (k'*2*pi).^2);
disp(sprintf('Mean energy dissipation rate = %d',epsilon));
