%
%########################################################################
%#  plot of DNS helicity spectrum output file
%########################################################################
%

%name = 'helicity_data/check256_hq_0030.0000';
namedir = '/home2/skurien/';
%name = 'helicity_data/sc1024_data/sc1024A0002.0000';
%mu = 0.35e-4;

name = 'dns/src/sk_helical0000.0000';
mu = 5.0e-4;


% plot all the spectrum:
movie=1;


spec_r_save=[];
spec_r_save_fac3=[];

fid=fopen([namedir,name,'.hspec'],'r','l');


[time,count]=fread(fid,1,'float64');
while (time >=.0 & time <=9999.3)
if (count==0) 
   disp('error reading hspec file')
end
time
  n_r=fread(fid,1,'float64');
  hspec_n=fread(fid,n_r,'float64');
  hspec_p = fread(fid,n_r,'float64');
scaling = 1;                       %2pi for Takeshi data, 1 otherwise

  hspec_n = hspec_n*scaling;
hspec_p = hspec_p*scaling;


figure(20);
loglog(abs(hspec_n),'r');  
hold on;
loglog(hspec_p,'b');
loglog(abs(hspec_p+hspec_n),'k');
loglog(abs(hspec_p-hspec_n),'k.-');
grid
%hold off

k = [1:n_r];
figure(21)
loglog(abs(hspec_n).*k'.^(4/3),'r');  
hold on;
loglog(hspec_p.*k'.^(4/3),'b');
loglog(abs(hspec_p-hspec_n).*k'.^(4/3),'k.-');
loglog(abs(hspec_p+hspec_n).*k'.^(4/3),'k');
grid on;
%hold off

figure(22)
loglog(abs(hspec_n).*k'.^(5/3),'r');  
hold on;
loglog(hspec_p.*k'.^(5/3),'b');
loglog(abs(hspec_p-hspec_n).*k'.^(5/3),'k.-');
loglog(abs(hspec_p+hspec_n).*k'.^(5/3),'k');
grid on;
%hold off


%compute mean helicity
H = sum(hspec_p+hspec_n)

%compute dissipation rate of helicity
h = sum((hspec_p+hspec_n).*k'.^2)*2*mu
[time,count]=fread(fid,1,'float64');
end
