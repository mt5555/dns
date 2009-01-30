%
%########################################################################
%#  plot of Craya-Herring projection spectra from *.spec_ch
%########################################################################


namedir ='/nh/u/skurien/projects/helicity/extract_helicity/';
name = 'skhel512a0007.0000';
mu = 1e-4;


% plot all the spectrum:
movie=0;


spec_r_save=[];
spec_r_save_fac3=[];

% note: endianopen() will notwork with .spec and  .hspec files
%because first number is not necessaryly an integer 
%pname = [strrep(name,'_','-'),'.hspec'];
fid=fopen([namedir,name,'.spec_ch'],'r','l');


spec_tot = [];
spec_wave = [];
spec_vort = [];
spec_kh0 = [];

[time,count]=fread(fid,1,'float64');
j = 0;
while (time >=.0 & time <= 9999.3)
if (count==0) 
   disp('error reading spec_ch file')
end

time
  n_r=fread(fid,1,'float64');
  spec_tot=fread(fid,n_r,'float64');
  spec_vort=fread(fid,n_r,'float64');
  spec_wave=fread(fid,n_r,'float64');
  spec_kh0=fread(fid,n_r,'float64');

%if (j == 0) 
%if (time < 1)
%     spec_ave = 0*(hspec_n + hspec_p); % initialize hspec_ave
%end
%if (time>=1)
%if (time >= 1)
%j=j+1
%if (length(hspec_ave)==0) 
%  hspec_ave = hspec_n + hspec_p;
%else
%  hspec_ave = hspec_ave + hspec_n + hspec_p;
%end
%end

k = [0:n_r-1];

if (movie)
figure(1); % +, - and total helicity spectra
loglog(k,spec_tot,'r'); hold on;
loglog(k,spec_vort,'b'); hold on;
loglog(k,spec_wave,'k'); hold on;
loglog(k,spec_kh0,'c'); 
%axis([1 100 1e-2 10]);
grid
legend('total','vortical','wave','kh=0')
%hold off

end % end movie loop

[time,count]=fread(fid,1,'float64');

end

