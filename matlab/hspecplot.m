%
%########################################################################
%#  plot of DNS helicity spectrum output file
%########################################################################
% make sure viscosity is in Marks units. All the .hspec files (even those from 
% Takeshi's data are in Marks units already

%name = 'helicity_data/helicity_spc/check256_hq_0000.0000';
%namedir = '/home2/skurien/';
%mu = 6.0d-3;                          % viscosity in Takeshi's units
%scaling = 2*pi;                       %2pi for Takeshi data, 1 otherwise;
%mu = mu/scaling^2;

%name = 'helicity_data/sc1024_data/sc1024A0002.0000';
%mu = 0.35e-4;

%name = 'dns/src/hel128_hpi4/hel128_hpi4_0000.0000';
%mu = 5.0e-4;

%name = 'dns/src/hel128_hpi2/hel128_hpi2_0000.0000';
%mu = 5.0e-4;

namedir = '/home2/skurien/helicity_data/helical_forced/';
name = 'hel256_hpi2/hel256_hpi2_all';
mu = 2e-4;

%name = 'dns/src/hel128_hpi2/hel128_hpi2_0000.0000';
%mu = 5e-4;


% plot all the spectrum:
movie=0;


spec_r_save=[];
spec_r_save_fac3=[];

fid=fopen([namedir,name,'.hspec'],'r','l');


hspec_n = [];
hspec_p = [];
hspec_ave = hspec_n + hspec_p;

[time,count]=fread(fid,1,'float64');
j = 0;
while (time >=.0 & time <= 9999.3)
j=j+1
if (count==0) 
   disp('error reading hspec file')
end
time
  n_r=fread(fid,1,'float64');
  hspec_n=fread(fid,n_r,'float64');
hspec_p = fread(fid,n_r,'float64');

if (j == 1) 
     hspec_ave = hspec_n + hspec_p;
else
hspec_ave = hspec_ave + hspec_n + hspec_p;
end

k = [0:n_r-1];

if (movie)
figure(20); % +, - and total helicity spectra
loglog(k,abs(hspec_n),'r'); hold on;
loglog(k,hspec_p,'b'); hold on;
loglog(k,abs(hspec_p+hspec_n),'k'); hold on;
%axis([1 100 1e-2 10]);
grid
legend('|H_-(k)|','H_+(k)','H(k)=H_+(k) + H_-(k)')
%hold off


figure(21); % 4/3 normalized +, - and total helicity spectra
loglog(k,abs(hspec_n).*k'.^(4/3),'r'); hold on; 
loglog(k,hspec_p.*k'.^(4/3),'b');
loglog(k,abs(hspec_p-hspec_n).*k'.^(4/3),'k-');hold on;
%loglog(k,abs(hspec_p+hspec_n).*k'.^(4/3),'k.-');hold on;
%axis([1 100 1e-2 10]);hold on;
grid on;
title('4/3 compensated helicity spectra');
legend('|H_-(k)| k^{4/3}','H_+(k) k^{4/3}','H(k) k^{4/3}')
%hold off

figure(22); % 5/3 normalized +, - and total helicity spectra
loglog(k,abs(hspec_n).*k'.^(5/3),'r');hold on;
loglog(k,hspec_p.*k'.^(5/3),'b');hold on;
loglog(k,abs(hspec_p-hspec_n).*k'.^(5/3),'k-');hold on;
%loglog(k,abs(hspec_p+hspec_n).*k'.^(5/3),'k.-'); hold on;
%axis([1 100 1e-2 10]);
grid on;
title('5/3 compensated helicity spectra');
legend('|H_-(k)| k^{5/3}','H_+(k) k^{5/3}','H(k) k^{5/3}')
%hold off

end % end movie loop

% compute total helicity for snapshot -- in Mark's units
% multiply by 2*pi^2 to get Takeshi's units
H = sum(hspec_p+hspec_n);
disp(sprintf('Total helicity = %d',H))

% compute dissipation rate of helicity -- in Mark's units 
% multiply dissipation rate by 2*pi to get Takeshi's units

h = -2*mu*sum((hspec_p+hspec_n).*(k'*2*pi).^2);
disp(sprintf('Helicity dissipation rate = %d',h));

figure(23)
plot(time,H,'o',time,h,'x');hold on;
legend('Total Helicity','Mean dissipation rate')
[time,count]=fread(fid,1,'float64');
end

hspec_ave = hspec_ave/(j);

figure(24)
loglog(k,abs(hspec_ave),'x'); hold on;
title('Average helicity spectrum')


figure(25)
loglog(k,abs(hspec_ave).*k'.^(4/3),'k');hold on;
title('Average 4/3 compensated helicity spectrum');

figure(26)
loglog(k,abs(hspec_ave).*k'.^(5/3),'b');hold on;
title('Average 5/3 compensated helicity spectrum');

% compute MEAN helicity in time -- in Mark's units
% multiply by 2*pi^2 to get Takeshi's units
H = sum(hspec_ave);
disp(sprintf('Total helicity = %d',H))

% compute dissipation rate of helicity -- in Mark's units 
% multiply dissipation rate by 2*pi to get Takeshi's units

h = -2*mu*sum((hspec_ave).*(k'*2*pi).^2);
disp(sprintf('Helicity dissipation rate = %d',h));


%end
