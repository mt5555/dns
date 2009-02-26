%
%########################################################################
%#  plot of Craya-Herring projection spectra from *.spec_ch
%########################################################################


namedir ='~kurien/INCITE_runs/SW02_tests/';
name = 'bous128_Ro21Fr0.21_all';


%namedir ='~kurien/dns/src/';
%name = 'temp_all';
%mu = 1e-4;


% plot all the spectrum:
movie=1;


spec_r_save=[];
spec_r_save_fac3=[];

fid=fopen([namedir,name,'.spec_ch'],'r','l');


spec_tot = [];
spec_Q_tot = [];
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
  spec_Q_tot=fread(fid,n_r,'float64');
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
%loglog(k,spec_Q_tot,'bo'); hold on;
%loglog(k,spec_tot./spec_Q_tot,'ko');hold on;
loglog(k,spec_vort,'b'); hold on;
loglog(k,spec_wave,'k'); hold on;
loglog(k,spec_kh0,'c');hold on; 
axis([1 1000 1e-6 1]);
grid
legend('total','vortical','wave','kh=0')
hold off
pause
end % end movie loop

[time,count]=fread(fid,1,'float64');

end

