%
%########################################################################
%#  plot of DNS helicity spectrum output file
%########################################################################
%

name = 'helicity_data/helicity_spc/check256_hq_0000.0000';
namedir = '/home2/skurien/';
%name = 'helicity_data/sc1024_data/sc1024A0002.0000';
%mu = 0.35e-4;

%name = 'dns/src/hel128_hpi4/hel128_hpi4_0000.0000';
%mu = 5.0e-4;

%name = 'dns/src/hel128_hpi2/hel128_hpi2_0000.0000';
%mu = 5.0e-4;

%name = 'dns/src/hel256_hpi2_0000.0000';
%mu = 2e-4;

name = 'dns/src/hel128_hpi2/hel128_hpi2_0000.0000';
mu = 5e-4;

% plot all the spectrum:
movie=1;


spec_r_save=[];
spec_r_save_fac3=[];

fid=fopen([namedir,name,'.hspec'],'r','l');


[time,count]=fread(fid,1,'float64');
j = 0;
while (time >=.0 & time <=9999.3)
if (count==0) 
   disp('error reading hspec file')
end
time
  n_r=fread(fid,1,'float64');
  hspec_n=fread(fid,n_r,'float64');
hspec_p = fread(fid,n_r,'float64');
scaling = 1;                       %2pi for Takeshi data, 1 otherwise;

hspec_n = hspec_n*scaling;
hspec_p = hspec_p*scaling;
k = [1:n_r];
%if (time == 5.0) then
figure(20);
%loglog(abs(hspec_n),'r'); hold on;
%hold on;
%loglog(hspec_p,'b'); hold on;
loglog(abs(hspec_p+hspec_n),'k'); hold on;
loglog(abs(hspec_p-hspec_n),'k.-'); hold on;
%axis([1 100 1e-2 10]);
grid
%hold off


figure(21);
%loglog(abs(hspec_n).*k'.^(4/3),'r'); hold on; 
%hold on;
%loglog(hspec_p.*k'.^(4/3),'b');
loglog(abs(hspec_p-hspec_n).*k'.^(4/3),'k.-');hold on;
loglog(abs(hspec_p+hspec_n).*k'.^(4/3),'k');hold on;
%axis([1 100 1e-2 10]);hold on;
grid on;
%hold off

figure(22); 
%loglog(abs(hspec_n).*k'.^(5/3),'r');hold on;
%hold on;
%loglog(hspec_p.*k'.^(5/3),'b');hold on;
loglog(abs(hspec_p-hspec_n).*k'.^(5/3),'k.-');hold on;
loglog(abs(hspec_p+hspec_n).*k'.^(5/3),'k'); hold on;
%axis([1 100 1e-2 10]);
grid on;
%hold off



%compute total helicity
H = sum(hspec_p+hspec_n)

%compute dissipation rate of helicity
h = -2*mu*sum((hspec_p+hspec_n).*(k').^2)
figure(23)
plot(time,H,'o',time,h,'x');hold on;
[time,count]=fread(fid,1,'float64');
j = j+1
end
%end
