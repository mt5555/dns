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

%name = 'skhel512a0009.2781';
%namedir = '/home/mt/data/skhel/';
namedir = '/nh/nest/u/skurien/projects/helicity_data/helical_forced/hel256_hpi2/';
name = 'hel256_hpi2_all';
mu = 2e-4;

%namedir = '/home2/skurien/helicity_data/helical_forced/hel128_h3pi8/';
%name = 'hel128_h3pi8_0000.0000';
%mu = 5e-4;

%namedir = '/home2/skurien/helicity_data/helical_forced/hel256_h0/';
%name = 'hel256_h0';
%mu = 2e-4;

namedir = '/nh/nest/u/skurien/projects/helicity_data/helical_forced/hel512_hpi2/diag/';
name = 'skhel512_hpi2';
mu = 1e-4;


%namedir = '/nh/nest/u/skurien/projects/helicity_data/helical_forced/hel512_hpi2/diag/';
%name = 'sc1024A';
%mu = 3.5e-5;



% plot all the spectrum:
movie=0;


spec_r_save=[];
spec_r_save_fac3=[];

% note: endianopen() will notwork with .spec and  .hspec files
%because first number is not necessaryly an integer 
pname = [strrep(name,'_','-'),'.hspec'];
fid=fopen([namedir,name,'.hspec'],'r','l');


hspec_n = [];
hspec_p = [];
hspec_ave = hspec_n + hspec_p;

[time,count]=fread(fid,1,'float64');
j = 0;
while (time >=.0 & time <= 9999.3)
if (count==0) 
   disp('error reading hspec file')
end
time
  n_r=fread(fid,1,'float64');
  hspec_n=fread(fid,n_r,'float64');
hspec_p = fread(fid,n_r,'float64');

%if (j == 0) 
if (time < 1)
     hspec_ave = 0*(hspec_n + hspec_p); % initialize hspec_ave
end
%if (time>=1)
if (time >= 1)
j=j+1
if (length(hspec_ave)==0) 
  hspec_ave = hspec_n + hspec_p;
else
  hspec_ave = hspec_ave + hspec_n + hspec_p;
end
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

%figure(23)
%plot(time,H,'o',time,h,'x');hold on;
%legend('Total Helicity','Mean dissipation rate')
%xlabel('time');
%ylabel(pname);
[time,count]=fread(fid,1,'float64');
end

hspec_ave = hspec_ave/(j);

figure(24)
loglog(k,abs(hspec_ave),'b'); hold on;
title('Average helicity spectrum')
xlabel('k')
ylabel(pname);

figure(25)
semilogx(k,abs(hspec_ave).*k'.^(5/3),'k');hold on;
semilogx(k,abs(hspec_ave).*k'.^(4/3),'b--');hold on;
set(gca,'fontsize',16);
%title('Average compensated helicity spectrum');
legend('H(k) k^{5/3}','H(k) k^{4/3}');
set(gca,'fontsize',18);	      
xlabel('k')
%ylabel(pname);


% compute MEAN helicity in time -- in Mark's units
% multiply by 2*pi^2 to get Takeshi's units
H = sum(hspec_ave);
disp(sprintf('Mean helicity = %d',H))

% compute dissipation rate of helicity -- in Mark's units 
% multiply dissipation rate by 2*pi to get Takeshi's units

h = -2*mu*sum((hspec_ave).*(k'*2*pi).^2);
disp(sprintf('Mean helicity dissipation rate = %d',h));


%end
