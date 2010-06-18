clear

forpaper = 0; % generate plots for paper
totals = 0; %general plots about total energies
scales = 1; %general plots showing evolution of correlation lengths
Lz = 1; %default


namedir = '~/projects/INCITE_runs/Intrepid/qg/';
name = 'n640_bous3000_all';
kf=4; Lz = 1;

%namedir = '~/projects/INCITE_runs/Intrepid/qg/n256_fatshellforce/';
%name = 'n256_Ro0.01_2_all';
%kf=16; Lz = 1;

namedir = '~/projects/INCITE_runs/Intrepid/lowaspect_bous/n2048_d0.25_Ro0.05_nodamp/';
name = 'n2048_d0.25_Ro0.05_all';
kf=4; Lz = 0.25;

%namedir = '~/projects/INCITE_runs/Intrepid/lowaspect_bous/n2048_d0.25_Ro0.01_nodamp/';
%name = 'n2048_d0.25_Ro0.01_all';
%kf=4; Lz = 0.25;

%namedir = '~/projects/INCITE_runs/Intrepid/lowaspect_bous/n2048_d0.25_Ro0.005_nodamp/';
%name = 'n2048_d0.25_Ro0.005_all';
%kf=4; Lz = 0.25;

namedir = '~/projects/INCITE_runs/Intrepid/lowaspect_bous/n2048_d0.25_Ro0.002_nodamp/';
name = 'n2048_d0.25_Ro0.002_all';
kf=4; Lz = 0.25;

asciprint = 0 % if == 1 print out the data to asci files

fid=fopen([namedir,name,'.spec_ch2d'],'r','b');  %use 'b' for 
if (fid<0) 
  'Error opening file',[namedir,name]
  return
end


time=0;
j = 0; %(count for time average)
while (time < 35)
  [time,count]=fread(fid,1,'float64');
  if (count ~= 1) 
    disp('EOF reached.  stopping')
    break 
  end
  if (time<0 | time>1000) break; end;
  numvar=fread(fid,1,'float64');
  numkh=fread(fid,1,'float64');
  numkz=fread(fid,1,'float64');
  spec2d_uwave=zeros([numkz,numkh]);  
  spec2d_vwave=zeros([numkz,numkh]);
  spec2d_uvort=zeros([numkz,numkh]);
  spec2d_vvort = zeros([numkz,numkh]);
  
  disp(sprintf('time=%f  kz=%f  kh=%f',time,numkz,numkh));

  for ivar=1:numvar
      for kz=1:numkz
        [s,count] = fread(fid,numkh,'float64') ;
     %   kz
     %   max(s)      
    if (count~=numkh)
      disp('Error: error reading file')
      count
      size(s)
    end
    if (ivar == 1)
        spec2d_uwave(kz,:) = s';
    elseif (ivar==2)
        spec2d_vwave(kz,:) = s';
    elseif (ivar==3)
        spec2d_uvort(kz,:) = s';
    elseif (ivar ==4)
        spec2d_vvort(kz,:) = s';
    end
    end
  end

  % kz-averaged 2d spectra (remaining array is function of kh)
  specuwave_h = sum(spec2d_uwave,1);
  specvwave_h = sum(spec2d_vwave,1);
  specuvort_h = sum(spec2d_uvort,1);
  specvvort_h = sum(spec2d_vvort,1);

  %kh-averaged 1d spectra (remaining array is function of kz)
  specuwave_z = sum(spec2d_uwave,2);
  specvwave_z = sum(spec2d_vwave,2);
  specuvort_z = sum(spec2d_uvort,2);
  specvvort_z = sum(spec2d_vvort,2);

  expo = 0 ; %0 for no compensation of spectra

  kz = (1:numkz)-1;
  kh = (1:numkh)-1;
  kzvals = [1,2,3,4,5,11,21,31,41,51];
  khvals = [1,2,3,4,5,11,21,31,41,51];

  %plot movie of spectra evolving in time
if(0)
  %plot E_h(kh,kz) and P(kh,kz) as a function of kh for various kz
  figure(1);
  set(gca,'fontsize',16);
  subplot(2,1,1);hold off;
  subplot(2,1,2); hold off;
  for i = 1:length(kzvals)
  subplot(2,1,1); 
  loglog(kh,(spec2d_uwave(kzvals(i),:)'+spec2d_vwave(kzvals(i),:)').*(kh'.^expo)); hold on;%pause
  ylabel('E^{\pm}_h(k_h,k_z)')
  subplot(2,1,2); 
  loglog(kh,spec2d_uvort(kzvals(i),:)'+spec2d_vvort(kzvals(i),:).*(kh'.^expo));hold on;%pause
  ylabel('E^0_h(k_h,k_z)')
  figure(1); xlabel('k_h')
  end  
 

  %plot E_h(kh,kz) and P(kh,kz) as a function of kz for various kh
  figure(2);
  set(gca,'fontsize',16);
  subplot(2,1,1); hold off;
  subplot(2,1,2); hold off;
  for i = 1:length(khvals)
  subplot(2,1,1)
  loglog(kz,(spec2d_uwave(:,khvals(i))+spec2d_vwave(:,khvals(i))).*(kz'.^(expo)));hold on;%pause
  ylabel('E^{\pm}_h(k_h,k_z)')
  subplot(2,1,2)
  loglog(kz,spec2d_uvort(:,khvals(i))+spec2d_vvort(:,khvals(i)).*(kz'.^(expo)));hold on;%pause
  ylabel('E^0_h(k_h,k_z)')
  figure(2); xlabel('k_z')
  end
  

end

if (1)
  %plot the kz-averaged spectra
  figure(3)
  %loglog53(numkh,specuwave_h'+specvwave_h','E^{\pm}_h(kh)',2.0,6);%hold on;
  loglog(kh,specuwave_h'+specvwave_h');
  ylabel('E^{\pm}_h(kh)');
  xlabel('k_h');
  figure(5)
  loglog(kh,specuvort_h'+specvvort_h');
  ylabel('E^0_h(kh)');%hold on;
  xlabel('k_h');
  
  %plot the kh-averaged spectra
  figure(6)
  loglog(kz/Lz,(specuwave_z+specvwave_z)*Lz);
  ylabel('E^{\pm}_h(kz)');
  xlabel('k_z')
  figure(7)
  loglog(kz/Lz,(specuvort_z + specvvort_z)*Lz);
  ylabel('E^0_h(kz)');
  xlabel('k_z')
end
  
  if (asciprint == 1)
    tstamp=10000+time;
    tstamp=sprintf('%.4f',tstamp);
    tstamp=tstamp(2:10);
    aname=sprintf('%s%s-%s.gnu',namedir,name,tstamp);
    fidd=fopen(aname,'w');
    if (fidd<0) 
      'Error opening file',aname
      return
    end
    for k=1:numkh
      fprintf(fidd,'%25.15e  %25.15e \n',k,spec2d(1,k));
    end 
    fprintf(fidd,'\n');
    fclose(fidd);
  
    aname=sprintf('%s%s-%s.ignu',namedir,name,tstamp);
    fidd=fopen(aname,'w');
    if (fidd<0) 
      'Error opening file',aname
      return
    end
    for k=1:numkh
      fprintf(fidd,'%25.15e  %25.15e \n',k,spec(k));
    end 
    fprintf(fidd,'\n');
    fclose(fidd);
  
  end
% disp('pause') 
% pause
 

% compute vertical and horizontal length scales of wave velocity and plot
if (scales == 1)
n = 1;
Ehkh_wave = sum(specuwave_h + specvwave_h);
Ehkz_wave = sum(specuwave_z + specvwave_z);
khEh = sum((2*pi*kh).^(1/n).*(specuwave_h + specvwave_h));
%khEh = sum((kh).^(1/n).*(specuwave_h + specvwave_h));
kzEh = sum((2*pi*kz/Lz).^(1/n)*(specuwave_z + specvwave_z));
%kzEh = sum((kz/Lz).^(1/n)*(specuwave_z + specvwave_z));
kH = (kzEh/Ehkz_wave)^n
kL = (khEh/Ehkh_wave)^n
H = 2*pi/kH;
L = 2*pi/kL;
delta = H/L
figure(20)
plot(time, H, 'x'); hold on;
set(gca,'fontsize',16);
xlabel('time');
ylabel('internal vertical scale H (wave)');
figure(21)
plot(time, L, 'o'); hold on;
set(gca,'fontsize',16);
xlabel('time');
ylabel('internal horizontal scale L (wave)');
figure(22)
plot(time, delta, '*'); hold on;
set(gca,'fontsize',16);
xlabel('time');
ylabel('internal scale aspect ratio H/L (wave)');
figure(23)
plot(time,Lz/H, 'd'); hold on;
set(gca,'fontsize',16);
xlabel('time');
ylabel('\delta/H');
end

% compute vertical and horizontal length scales of vortical velocity and plot
if (scales == 1)
n = 1;
Ehkh_vort = sum(specuvort_h + specvvort_h);
Ehkz_vort = sum(specuvort_z + specvvort_z);
khEh = sum((2*pi*kh).^(1/n).*(specuvort_h + specvvort_h));
kzEh = sum((2*pi*kz/Lz).^(1/n)*(specuvort_z + specvvort_z));
kH = (kzEh/Ehkz_vort)^n;
kL = (khEh/Ehkh_vort)^n;
H = 2*pi/kH;
L = 2*pi/kL;
delta = H/L
figure(24)
plot(time, H, 'x'); hold on;
set(gca,'fontsize',16);
xlabel('time');
ylabel('internal vertical scale H (vortical)');
figure(25)
plot(time, L, 'o'); hold on;
set(gca,'fontsize',16);
xlabel('time');
ylabel('internal horizontal scale L (vortical)');
figure(26)
plot(time, delta, '*'); hold on;
set(gca,'fontsize',16);
xlabel('time');
ylabel('internal scale aspect ratio, H/L (vortical)');
figure(27)
plot(time,Lz/H, 'd'); hold on;
set(gca,'fontsize',16);
xlabel('time');
ylabel('\delta/H (vortical)');
end



 %time average of spectra
 if(0)
 if (time > 0 & time <= 5)
    if (j==0) 
        spec2d_Ehwave_ave = spec2d_u + spec2d_v;
        spec2d_Ev_ave = spec2d_w;
        specEhkz_ave = specuz + specvz; %avg over kh
        specEvkz_ave = specwz + specwz;
        specPkz_ave = spectz;           %avg over kh
        specEhkh_ave = specuh + specvh; %avg over kz
        specEvkh_ave = specwh;
        specPkh_ave = specth;           %avg over kz        
        j = j+1;
    else
        spec2d_t_ave = spec2d_t_ave + spec2d_t;
        spec2d_Eh_ave = spec2d_Eh_ave + spec2d_u + spec2d_v;
        specEhkz_ave = specEhkz_ave + specuz + specvz;
        spec2d_Ev_ave = spec2d_Ev_ave + spec2d_w;
        specPkz_ave = specPkz_ave+ spectz;
        specEhkh_ave = specEhkh_ave + specuh + specvh;
        specEvkh_ave = specEvkh_ave + specwh;
        specPkh_ave = specPkh_ave + specth;
        
        j = j + 1;
    end
 end
 end 
end
if(0)
spec2d_t_ave = spec2d_t_ave/j;
spec2d_Eh_ave = spec2d_Eh_ave/j;
specEhkz_ave = specEhkz_ave/j;
specEvkz_ave = specEvkz_ave/j;
specPkz_ave = specPkz_ave/j;
specEhkh_ave = specEhkh_ave/j;
specEvkh_ave = specEvkh_ave/j;
specPkh_ave = specPkh_ave/j;

%plot time avg spectra

figure(10);hold off;
subplot(2,1,1);hold off;
subplot(2,1,2); hold off;
for i = 1:length(khvals)     
subplot(2,1,1) 
axis([1 640/3 1e-10 1])
set(gca,'fontsize',17);
loglog(kz,spec2d_Eh_ave(:,khvals(i)).*((kz'.^expo)));hold on;%pause;
ylabel('E_h(k_h,k_z)')
subplot(2,1,2)
axis([1 640/3 1e-10 1])
set(gca,'fontsize',17);
loglog(kz,spec2d_t_ave(:,khvals(i)).*((kz'.^expo)));hold on;%pause;
ylabel('P(k_h,k_z)');
end
figure(10);subplot(2,1,2)
xlabel('k_z');
   
figure(11);hold off; 
subplot(2,1,1);set(gca,'fontsize',17); hold off;
subplot(2,1,2);set(gca,'fontsize',17); hold off;
for i = 1:length(kzvals) 
subplot(2,1,1)
loglog(kh,spec2d_Eh_ave(kzvals(i),:).*((kh.^expo)));hold on;%pause
ylabel('E_h(k_h,k_z)')
subplot(2,1,2)
loglog(kh,spec2d_t_ave(kzvals(i),:).*((kh.^expo)));hold on;%pause
ylabel('P(k_h,k_z)');
end
figure(11);
xlabel('k_h')

figure(12);hold off; 
subplot(2,1,1);hold off;
loglog(kh,specEhkh_ave.*((kh.^expo)));hold on;%pause
set(gca,'fontsize',17)
ylabel('E_h(k_h)');
subplot(2,1,2);hold off;
loglog(kh,specPkh_ave.*((kh.^expo)));hold on;%pause
set(gca,'fontsize',17)
ylabel('P(k_h)');
xlabel('k_h');
figure(12);subplot(2,1,1);title('time ave');

figure(13); hold off;
subplot(2,1,1);hold off;
loglog(kz,specEhkz_ave.*((kz'.^expo)));hold on;%pause
set(gca,'fontsize',17)
ylabel('E_h(k_z)');
subplot(2,1,2);hold off;
loglog(kz,specPkz_ave.*((kz'.^expo)));hold on;%pause
set(gca,'fontsize',17)
ylabel('P(k_z)');
xlabel('k_z')
figure(13);subplot(2,1,1); title('Time ave')

figure(14); hold off;
subplot(2,1,1);hold off;
loglog(kz,specEvkz_ave.*((kz'.^expo)));hold on;%pause
set(gca,'fontsize',17)
ylabel('E_v(k_z)');
xlabel('k_z');
title('Time ave');
subplot(2,1,2);hold off;
loglog(kh,specEvkh_ave.*((kh.^expo)));hold on;%pause
set(gca,'fontsize',17);
ylabel('E_v(k_h)');
xlabel('k_h');


%normalize x-axis by ratio kz or kh as necessary
figure(20)
for i = 1:length(khvals)
loglog(kz/khvals(i),spec2d_t_ave(:,khvals(i)));hold on;%pause
xlabel('$kz/kh$')
ylabel('$P(k_h,k_z)$')
title('P(kh,kz) vs. kz/kh')
end

figure(21)
for i = 1:length(kzvals)
loglog(kh/kzvals(i),spec2d_Eh_ave(kzvals(i),:));hold on;%pause
xlabel('$kh/kz$')
ylabel('$E_h(k_h,k_z)$')
title('Eh(kh,kz) vs. kh/kz')

end
end

if (forpaper == 1) 
 
figure(30);hold off;
for i = 1:length(khvals)     
axis([1 640/3 1e-10 10^5])
set(gca,'fontsize',17);
loglog(kz,10e4*spec2d_Eh_ave(:,khvals(i)).*((kz'.^expo)));hold on;%pause;
loglog(kz,spec2d_t_ave(:,khvals(i)).*((kz'.^expo)),'r--');hold on;%pause;
legend('E_h(k_h,k_z) \times 10^4','P(k_h,k_z)');
x = 5:40;
y = x.^(-5);
plot(x,y,'b');
xx = 10:80;
yy = 10^5*xx.^(-3);
plot(xx,yy,'r');
end
xlabel('k_z');
   
figure(31);hold off;
for i = 1:length(kzvals) 
axis([1 640/3 1e-15 10^5])
set(gca,'fontsize',17);
loglog(kh,10e4*spec2d_Eh_ave(kzvals(i),:).*((kh.^expo)));hold on;%pause
loglog(kh,spec2d_t_ave(kzvals(i),:).*((kh.^expo)),'r--');hold on;%pause
legend('E_h(k_h,k_z) \times 10^4','P(k_h,k_z)');
x = 5:40;
y = 10^5*x.^(-4.5);
plot(x,y,'b');
xx = 10:80;
yy = xx.^(-2);
plot(xx,yy,'r');
end
xlabel('k_h')

figure(32);hold off; 
axis([1 640/3 1e-20 10])
set(gca,'fontsize',17)
loglog(kh,specEhkh_ave.*((kh.^expo)),'Linewidth',2);hold on;%pause
loglog(kh,specPkh_ave.*((kh.^expo)),'r--','Linewidth',2);hold on;%pause
legend('E_h(k_h)','P(k_h)');
x = 5:40;
y = x.^(-1);
plot(x,y,'k');
xx = 10:80;
yy = xx.^(-4/3);
plot(xx,yy,'r');
xlabel('k_h');


figure(33); hold off;
axis([1 640/3 1e-10 10])
set(gca,'fontsize',17)
loglog(kz,specEhkz_ave.*((kz'.^expo)),'Linewidth',2);hold on;%pause
loglog(kz,specPkz_ave.*((kz'.^expo)),'r--','Linewidth',2);hold on;%pause
legend('E_h(k_z)','P(k_z)');
x = 5:40;
y = x.^(-1);
plot(x,y,'k');
xx = 10:80;
yy = xx.^(-4/3);
plot(xx,yy,'r');
xlabel('k_z')
end


if (totals == 1)

figure(42);hold off; 
axis([1 640/3 1e-20 10])
set(gca,'fontsize',17)
loglog(kh,(specEhkh_ave+specPkh_ave).*((kh.^expo)));hold on;%pause
legend('E_h(k_h)+P(k_h)');
x = 5:40;
y = x.^(-1);
plot(x,y,'k');
xx = 10:80;
yy = xx.^(-4/3);
plot(xx,yy,'r');
xlabel('k_h');


figure(43); hold off;
axis([1 640/3 1e-10 10])
set(gca,'fontsize',17)
loglog(kz,(specEhkz_ave + specPkz_ave).*((kz'.^expo)));hold on;%pause
legend('E_h(k_z)+P(k_z)');
x = 5:40;
y = x.^(-1);
plot(x,y,'k');
xx = 10:80;
yy = xx.^(-4/3);
plot(xx,yy,'r');
xlabel('k_z')  

end