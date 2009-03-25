%
%########################################################################
%#  plot of Craya-Herring projection spectra from *.spec_ch
%########################################################################

namedir ='~kurien/INCITE_runs/SW02_tests/bous128_Fr0.21/';
name = 'bous128_Fr0.21_all';

%namedir ='~kurien/INCITE_runs/SW02_tests/bous128_Ro21Fr0.21/';
%name = 'bous128_Ro21Fr0.21_all';

%namedir ='~kurien/INCITE_runs/SW02_tests/bous128_Ro2.1Fr0.21/';
%name = 'bous128_Ro2.1Fr0.21_all';


%namedir ='~kurien/INCITE_runs/RemSukSmi09_tests/lowres/';
%name = 'RSSlowres_all';

namedir ='~kurien/INCITE_runs/RemSukSmi09_tests/highres/';
name = 'RSShighres_all';

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


k = [0:n_r-1];

if (movie)
figure(1); % +, - and total and projected energy spectra
loglog(k,spec_tot,'r'); hold on;
%loglog(k,spec_Q_tot,'bo'); hold on;
%loglog(k,spec_tot./spec_Q_tot,'ko');hold on;
loglog(k,spec_vort,'b'); hold on;
loglog(k,spec_wave,'k'); hold on;
loglog(k,spec_kh0,'c');hold on; 
axis([1 1000 1e-6 1]);
grid
legend('total','vortical','wave')
hold off
pause

figure(2) ;
te = sum(spec_tot);
tvort = sum(spec_vort);
twave = sum(spec_wave);
tsk = (0.5*(2*pi*4)^2)^(1/3);
tls = (1*4^2)^(1/3);
esk = (0.5/2/pi/4)^(-2/3);
els = (1/4)^(-2/3);
plot(time*tsk/tls, te*esk/els,'kx'); hold on;
plot(time*tsk/tls, tvort*esk/els,'b.'); hold on;
plot(time*tsk/tls, twave*esk/els,'ro'); hold on;
legend('total','vortica','wave');

end % end movie loop

[time,count]=fread(fid,1,'float64');

end

