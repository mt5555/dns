%
%########################################################################
%#  plot of DNS 'tuned' helicity spectrum output file (*.hf_spec)
%########################################################################
% make sure viscosity is in Marks units. All the .hf_free files 
%(even those from Takeshi's data) are in Marks units already


namedir ='/netscratch/skurien/projects/helicity_data/helical_forced/hel256_hpi2/';
name = 'hel256_hpi2_0009.0000';
mu = 2e-4;

%namedir = '/home2/skurien/helicity_data/helical_forced/hel256_h0/';
%name = 'hel256_h0';
%mu = 2e-4;

%namedir = '/netscratch/skurien/projects/helicity_data/helical_forced/hel512_hpi2/diag/';
%name = 'skhel512_hpi2';
%mu = 1e-4;


%namedir = '/netscratch/skurien/projects/helicity_data/helical_forced/hel512_hpi2/diag/';
%name = 'sc1024A';
%mu = 3.5e-5;



% plot all the spectrum:
movie=0;


spec_r_save=[];
spec_r_save_fac3=[];

% note: endianopen() will notwork with .spec and  .hspec files
%because first number is not necessaryly an integer 
pname = [strrep(name,'_','-'),'.hf_spec'];
fid=fopen([namedir,name,'.hf_spec'],'r','l');

costta_spec = [];
hspec_n = [];
hspec_p = [];
hspec_ave = hspec_n + hspec_p;
espec = [];
kekspec = [];


[time,count]=fread(fid,1,'float64');
j = 0;
while (time >=.0 & time <= 9999.3)
if (count==0) 
   disp('error reading hf_spec file')
end
time
  n_r=fread(fid,1,'float64')
  hspec_n=fread(fid,n_r,'float64');
hspec_p = fread(fid,n_r,'float64');
espec = fread(fid,n_r,'float64');
kekspec = fread(fid,n_r,'float64');
costta_spec = fread(fid,n_r,'float64');
nbin = fread(fid,1,'int');
costta_pdf = zeros(n_r,nbin);
for i = 1:n_r
costta_pdf(i,:) = (fread(fid,nbin,'float64'))';
end

% bin values for pdfs of cos_tta
minct = -1
maxct = 1
ibin = (maxct - minct)/nbin
for i = 1:nbin
binvals(i)=minct+ibin*i;
end


%if (j == 0) 
if (time < 1)
     hspec_ave = 0*(hspec_n + hspec_p); % initialize hf_spec_ave
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
semilogx(k,abs(hspec_ave).*k'.^(5/3),'k','linewidth',[2]);hold on;
semilogx(k,abs(hspec_ave).*k'.^(4/3),'b-.','linewidth',[2]);hold on;
set(gca,'fontsize',16);
%title('Average compensated helicity spectrum');
legend('H(k) k^{5/3}','H(k) k^{4/3}');
set(gca,'fontsize',18);	      
xlabel('k')
%ylabel(pname);

figure(26); % tta spec
semilogx(k,abs(costta_spec),'r');hold on;
grid on;
title('helicity angle spectrum');
%hold off

figure(27)
loglog(k,0.5*espec,'r');hold on;
grid on;
title('energy spectrum');

figure(28)
for i = 1:n_r
plot(binvals, costta_pdf(i,:),'r')
hold on;
end


% compute MEAN helicity in time -- in Mark's units
% multiply by 2*pi^2 to get Takeshi's units
H = sum(hspec_ave);
disp(sprintf('Mean helicity = %d',H))

% compute dissipation rate of helicity -- in Mark's units 
% multiply dissipation rate by 2*pi to get Takeshi's units

h = -2*mu*sum((hspec_ave).*(k'*2*pi).^2);
disp(sprintf('Mean helicity dissipation rate = %d',h));

% compute MEAN energy in time -- in Mark's units
E= sum(espec);
disp(sprintf('Mean energy = %d',E))

%end
