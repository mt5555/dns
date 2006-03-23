%
%########################################################################
%#  plot of DNS 'tuned' helicity spectrum output file (*.hf_spec)
%########################################################################
% make sure viscosity is in Marks units. All the .hf_free files 
%(even those from Takeshi's data) are in Marks units already

%namedir ='/nh/u/skurien/projects/helicity/extract_helicity/';
%name = 'sc1024A0002.5000';
%mu = .35e-4;

%namedir ='/nh/u/skurien/projects/helicity/extract_helicity/';
%name = 'skhel512a_7to8.8';
%mu = 1e-4;

%namedir ='/nh/u/skurien/projects/helicity/extract_helicity/';
%name = 'skhel512a0007.0000';
%mu = 1e-4;

namedir ='/nh/u/skurien/projects/helicity/extract_helicity/';
name = 'skhel512a_uperp_0007.0000';
mu = 1e-4;

%namedir ='/nh/u/skurien/projects/helicity/extract_helicity/';
%name = 'skhel512a_upar_0007.0000';
%mu = 1e-4;

%namedir ='/nh/u/skurien/projects/helicity/extract_helicity/';
%name = 'skhel512a_del0.1_0007.0000.new.';
%mu = 1e-4;

%namedir ='/nh/u/skurien/projects/helicity/extract_helicity/';
%name = 'hel256_hpi2_0004.2000_nb1000';
%mu = 2e-4;

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
costta_spec_ave = [];
cosphi_spec = [];
cosphi_spec_ave = [];
hspec_n = [];
hspec_p = [];
hspec_ave = hspec_n + hspec_p;
espec = [];
spec_ave = [];
kekspec = [];
kekspec_ave = [];
e33 = [];
e33_ave = [];
ii2 = [];
ii2_ave = [];
rr2 = [];
rr2_ave = [];
i2i3 = [];
i2i3_ave = [];
r2i3 = [];
r2i3_ave = [];
par = [];
par_ave = [];
costta_pdf_ave = [];
cosphi_pdf_ave = [];
rhel_rms_ave = [];
costta_pdf = [];


[time,count]=fread(fid,1,'float64')
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
e33 = fread(fid,n_r,'float64'); %square component of I orthogonal to R
ii2 = fread(fid,n_r,'float64'); %square component of I along R
rr2 = fread(fid,n_r,'float64'); %square of R
i2i3 = fread(fid,n_r,'float64'); 
r2i3 = fread(fid,n_r,'float64');
par = fread(fid,n_r,'float64'); %|R\crossI|/(R^2 + I^2)
costta_specn = fread(fid,n_r,'float64');
costta_specp = fread(fid,n_r,'float64');
cosphi_specn = fread(fid,n_r,'float64');
cosphi_specp = fread(fid,n_r,'float64');
rhel_rms = fread(fid,n_r,'float64');
nbin = fread(fid,1,'float64'); %nbin is an integer stored as a float
%nbin=100;
costta_pdf = zeros(n_r,nbin);
cosphi_pdf = zeros(n_r,nbin); %the relative helicity pdf
for i = 1:n_r
costta_pdf(i,:) = (fread(fid,nbin,'float64'))';
end
for i = 1:n_r
cosphi_pdf(i,:) = (fread(fid,nbin,'float64'))';
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

if (length(spec_ave)==0)  
  spec_ave = espec;
  costta_specn_ave = costta_specn;
  costta_specp_ave = costta_specp;
   cosphi_specn_ave = cosphi_specn;
  cosphi_specp_ave = cosphi_specp;
    kekspec_ave = kekspec;
    e33_ave = e33;
    ii2_ave = ii2;
    rr2_ave = rr2;
    i2i3_ave = i2i3;
    r2i3_ave = r2i3;
    par_ave = par;
    rhel_rms_ave = rhel_rms;
    
    for i = 1:n_r
      costta_pdf_ave(i,:) = costta_pdf(i,:);
      cosphi_pdf_ave(i,:) = cosphi_pdf(i,:);
    end
    
    
else
  spec_ave = spec_ave + espec;
  costta_specn_ave = costta_specn_ave + costta_specn;
  costta_specp_ave = costta_specp_ave + costta_specp;  
  cosphi_specn_ave = cosphi_specn_ave + cosphi_specn;
  cosphi_specp_ave = cosphi_specp_ave + cosphi_specp;  
  kekspec_ave = kekspec_ave + kekspec;
  e33_ave = e33_ave + e33;
  ii2_ave = ii2_ave + ii2;
  rr2_ave = rr2_ave + rr2;
  i2i3_ave = i2i3_ave+i2i3;
  r2i3_ave = r2i3_ave+r2i3;
  par_ave = par_ave+ par;
  rhel_rms_ave = rhel_rms_ave + rhel_rms;
  for i = 1:n_r
    costta_pdf_ave(i,:) = costta_pdf_ave(i,:)+ costta_pdf(i,:);
    cosphi_pdf_ave(i,:) = cosphi_pdf_ave(i,:)+cosphi_pdf(i,:);
  end
  
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

% and for energy
E = sum(espec);
disp(sprintf('Total energy = %d',E))

e = mu*sum(espec.*(k'*2*pi).^2);
disp(sprintf('Energy dissipation rate = %d',e))


%figure(23)
%plot(time,H,'o',time,h,'x');hold on;
%legend('Total Helicity','Mean dissipation rate')
%xlabel('time');
%ylabel(pname);
[time,count]=fread(fid,1,'float64');
end

hspec_ave = hspec_ave/(j);
spec_ave = spec_ave/(j);
costta_specn_ave = costta_specn_ave/(j);
costta_specp_ave = costta_specp_ave/(j);
cosphi_specn_ave = cosphi_specn_ave/(j);
cosphi_specp_ave = cosphi_specp_ave/(j);
e33_ave = e33_ave/(j);
ii2_ave = ii2_ave/(j);
rr2_ave = rr2_ave/(j);
i2i3_ave = i2i3_ave/(j);
r2i3_ave = r2i3_ave/(j);
par_ave = par_ave/(j);
rhel_rms_ave = rhel_rms_ave/(j);
for i = 1:n_r
  costta_pdf_ave(i,:) = costta_pdf_ave(i,:)/j;
  cosphi_pdf_ave(i,:) = cosphi_pdf_ave(i,:)/j;
end


figure(24)
loglog(k,abs(hspec_ave),'b'); hold on;
title('Average helicity spectrum')
xlabel('k')
ylabel(pname);

figure(25)
semilogx(k,abs(hspec_ave).*k'.^(5/3),'k','linewidth',[2]);hold on;
semilogx(k,abs(hspec_ave).*k'.^(4/3),'b-.','linewidth',[2]);hold on;
set(gca,'fontsize',16);
title('Average compensated helicity spectra');
legend('H(k) k^{5/3}','H(k) k^{4/3}');
set(gca,'fontsize',18);	      
xlabel('k')
%ylabel(pname);

figure(26); % tta spec
semilogx(k,(costta_specn_ave),'r');hold on;
semilogx(k,(costta_specp_ave),'k');hold on;
grid on;
title('Spectrum of cosine of helicity angle');
%hold off

figure(27)
%semilogx(k,(cosphi_specn_ave),'r');hold on;
%semilogx(k,(cosphi_specp_ave),'k');hold on;
semilogx(k,(cosphi_specp_ave+cosphi_specn_ave),'b');hold on;
%semilogx(k,abs(hspec_ave)./kekspec,'m');hold on;
semilogx(k,rhel_rms_ave./(cosphi_specp_ave+cosphi_specn_ave), 'k');
grid on;
title('Spectrum of relative helicity');
%hold off


figure(28)
for i = 10:n_r-20
semilogy(binvals, costta_pdf_ave(i,:),'-');
title('PDF of Cosine of angle between RR and II for each wavenumber')
xlabel('cos(\theta)')
ylabel('P(cos\theta)')
hold on;
%pause
end

figure(35)
for i = 10:n_r-20
plot(binvals, cosphi_pdf_ave(i,:),'-');
hold on;
title('PDF of relative helicity for each wavenumber')
xlabel('Relative helicity')
%pause
end


figure(29)
semilogx(k,0.5*spec_ave.*k'.^(5/3),'k','linewidth',[2]);hold on;
semilogx(k,0.5*spec_ave.*k'.^(4/3),'b-.','linewidth',[2]);hold on;
%semilogx(k,0.5*spec_ave,'k','linewidth',[2]);hold on;
%semilogx(k,0.5*spec_ave,'b-.','linewidth',[2]);hold on;

set(gca,'fontsize',16);
title('Average compensated energy spectra');
legend('E(k) k^{5/3}','E(k) k^{4/3}');
set(gca,'fontsize',18);	      
xlabel('k')

figure(33)
semilogx(k,0.5*e33_ave.*k'.^(5/3),'k','linewidth',[2]);hold on;
semilogx(k,0.5*e33_ave.*k'.^(4/3),'b-.','linewidth',[2]);hold on;
set(gca,'fontsize',16);
title('Average compensated energy spectra');
legend('I_3^2(k) k^{5/3}','I_3^2(k) k^{4/3}');
set(gca,'fontsize',18);	      
xlabel('k')

figure(34)
semilogx(k,0.5*ii2_ave.*k'.^(5/3),'k','linewidth',[2]);hold on;
semilogx(k,0.5*ii2_ave.*k'.^(4/3),'b-.','linewidth',[2]);hold on;
set(gca,'fontsize',16);
title('Average compensated energy spectra');
legend('I_2^2(k) k^{5/3}','I_2^2(k) k^{4/3}');
set(gca,'fontsize',18);	      
xlabel('k')

figure(36)
semilogx(k,0.5*rr2_ave.*k'.^(5/3),'k','linewidth',[2]);hold on;
semilogx(k,0.5*rr2_ave.*k'.^(4/3),'b-.','linewidth',[2]);hold on;
set(gca,'fontsize',16);
title('Average compensated energy spectra');
legend('RR^2(k) k^{5/3}','RR^2(k) k^{4/3}');
set(gca,'fontsize',18);	      
xlabel('k')

figure(37)
semilogx(k,0.5*i2i3_ave.*k'.^(5/3),'k','linewidth',[2]);hold on;
semilogx(k,0.5*i2i3_ave.*k'.^(4/3),'b-.','linewidth',[2]);hold on;
set(gca,'fontsize',16);
title('Average compensated energy spectra');
legend('I_2I_3(k) k^{5/3}','I_2I_3(k) k^{4/3}');
set(gca,'fontsize',18);	      
xlabel('k')

figure(38)
semilogx(k,0.5*r2i3_ave.*k'.^(5/3),'k','linewidth',[2]);hold on;
semilogx(k,0.5*r2i3_ave.*k'.^(4/3),'b-.','linewidth',[2]);hold on;
set(gca,'fontsize',16);
title('Average compensated energy spectra');
legend('R_2I_3(k) k^{5/3}','R_2I_3(k) k^{4/3}');
set(gca,'fontsize',18);	      
xlabel('k')

figure(39)
loglog(k,0.5*par_ave,'k','linewidth',[2]);hold on;
set(gca,'fontsize',16);
title('Average compensated energy spectra');
legend('|R\timesI|/(R^2+I^2)');
set(gca,'fontsize',18);	      
xlabel('k')


% compute MEAN helicity in time -- in Mark's units
% multiply by 2*pi^2 to get Takeshi's units
H = sum(hspec_ave);
disp(sprintf('Mean helicity = %d',H))

% compute dissipation rate of helicity -- in Mark's units 
% multiply dissipation rate by 2*pi to get Takeshi's units

h = -2*mu*sum((hspec_ave).*(k'*2*pi).^2);
disp(sprintf('Mean helicity dissipation rate = %d',h));

% compute MEAN energy in time -- in Mark's units
E= sum(spec_ave);
disp(sprintf('Mean energy = %d',E))

%compute MEAN energy dissipation in time -- in Mark's units
e = mu*sum(spec_ave.*(k'*2*pi).^2);
disp(sprintf('Mean energy dissipation rate = %d',e));

%figure(30)
%semilogx(k,((0.5*espec - 0.2*(e*abs(h^(-1/3)).*(2*pi*k').^(-4/3))))./(1.8*e^(2/3).*(2*pi*k').^(-5/3)),'k-.','linewidth',[2]);hold on;
%semilogx(k,((0.5*espec - 1.8*(e*(2/3)).*(2*pi*k').^(-5/3)))./(0.2*e*abs(h^(-1/3)).*(2*pi*k').^(-4/3)),'b-.','linewidth',[2]);hold on;
%set(gca,'fontsize',14);
%title('');
%%legend('E(k)-0.4*\epsilon^{2/3}* k^{-5/3}');
%set(gca,'fontsize',14);	      
%xlabel('k')

%figure(31)
%semilogx(k,0.5*espec.*(e)^(-2/3).*(2*pi*k)'.^(5/3),'k','linewidth',[2]);hold on;
%semilogx(k,0.5*espec.*(e^(-1))*abs(h^(1/3)).*(2*pi*k)'.^(4/3),'b-.','linewidth',[2]);hold on;
%set(gca,'fontsize',16);
%title('Average compensated energy spectra');
%legend('E(k) k^{5/3}','E(k) k^{4/3}');
%set(gca,'fontsize',18);	      
%xlabel('k')


%end
