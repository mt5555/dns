%
%########################################################################
%#  plot of DNS helicity spectrum output file
%########################################################################
%

name = 'helicity_data/helicity_spc/check256_hq_cat30';
namedir = '/home2/skurien/';
mu = 6.0d-3;

%name = 'helicity_data/sc1024_data/sc1024A0002.0000';
%mu = 0.35e-4;

%name = 'dns/src/hel128_hpi4/hel128_hpi4_0000.0000';
%mu = 5.0e-4;

%name = 'dns/src/hel128_hpi2/hel128_hpi2_0000.0000';
%mu = 5.0e-4;

%namedir = '/home2/skurien/helicity_data/helical_forced/';
%name = 'hel256_hpi2/hel256_hpi2_cat';
%mu = 2e-4;

%name = 'dns/src/hel128_hpi2/hel128_hpi2_0000.0000';
%mu = 5e-4;


% plot all the spectrum:
movie=1;


spec_r_save=[];
spec_r_save_fac3=[];

fid=fopen([namedir,name,'.hspec'],'r','l');


hspec_n = [];
hspec_p = [];
hspec_ave = hspec_n + hspec_p;

[time,count]=fread(fid,1,'float64');
j = 0;
while (time >=.0 & time <=9999.3)
j=j+1
if (count==0) 
   disp('error reading hspec file')
end
time
  n_r=fread(fid,1,'float64');
  hspec_n=fread(fid,n_r,'float64');
hspec_p = fread(fid,n_r,'float64');
scaling = 1;                       %2pi for Takeshi data, 1 otherwise;

if (j == 1) 
     hspec_ave = hspec_n + hspec_p;
else
hspec_ave = hspec_ave + hspec_n + hspec_p;
end

hspec_n = hspec_n*scaling;
hspec_p = hspec_p*scaling;
k = [1:n_r];
%if (time == 5.0) then
figure(20);
loglog(abs(hspec_n),'r'); hold on;
loglog(hspec_p,'b'); hold on;
loglog(abs(hspec_p+hspec_n),'k'); hold on;
%loglog(abs(hspec_p-hspec_n),'k.-'); hold on;
%axis([1 100 1e-2 10]);
grid
legend('|H_-(k)|','H_+(k)','H(k)=H_+(k) + H_-(k)')
%hold off


figure(21);
loglog(abs(hspec_n).*k'.^(4/3),'r'); hold on; 
loglog(hspec_p.*k'.^(4/3),'b');
loglog(abs(hspec_p-hspec_n).*k'.^(4/3),'k-');hold on;
%loglog(abs(hspec_p+hspec_n).*k'.^(4/3),'k.-');hold on;
%axis([1 100 1e-2 10]);hold on;
grid on;
title('4/3 compensated helicity spectra');
legend('|H_-(k)| k^{4/3}','H_+(k) k^{4/3}','H(k) k^{4/3}')
%hold off

figure(22); 
loglog(abs(hspec_n).*k'.^(5/3),'r');hold on;
loglog(hspec_p.*k'.^(5/3),'b');hold on;
loglog(abs(hspec_p-hspec_n).*k'.^(5/3),'k-');hold on;
%loglog(abs(hspec_p+hspec_n).*k'.^(5/3),'k.-'); hold on;
%axis([1 100 1e-2 10]);
grid on;
title('5/3 compensated helicity spectra');
legend('|H_-(k)| k^{5/3}','H_+(k) k^{5/3}','H(k) k^{5/3}')
%hold off



%compute total helicity
H = sum(hspec_p+hspec_n);
disp(sprintf('Total helicity = %d',H))

%compute dissipation rate of helicity -- BUT THIS IS WRONG
h = -2*mu*sum((hspec_p+hspec_n).*(k').^2);
disp(sprintf('Mean helicity dissipation rate = %d',h));

figure(23)
plot(time,H,'o',time,h,'x');hold on;
legend('Total Helicity','Mean dissipation rate')
[time,count]=fread(fid,1,'float64')
%time,count
%j = j+1
end

hspec_ave = hspec_ave/(j-3);

figure(24)
loglog(abs(hspec_ave),'x'); hold on;

figure(25)
loglog(abs(hspec_ave).*k'.^(4/3),'k');hold on;

figure(26)
loglog(abs(hspec_ave).*k'.^(5/3),'b');hold on;


%end
