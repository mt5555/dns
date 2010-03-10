clear

forpaper = 1; % generate plots for paper
totals = 1; %general plots about total energies

epsilon=.41;
%CK=1.5*epsilon^(2/3);
%namedir = '/home/wingate/data1/Rotation/r16c/';

%namedir = '~/projects/pv/data_analysis/lowforc/low4/qg/qg64/sto_high_4/hyper_nu/bous100/';
%namedir = ['~/projects/pv/data_analysis/lowforc/low4/qg/qg64/' ...
%           'sto_high_4/hyper_nu/bous500/'];
%namedir = ['~/projects/pv/data_analysis/lowforc/low4/qg/qg64/' ...
%           'sto_high_4/hyper_nu/bous1000/'];
%namedir = '~/projects/pv/data_analysis/lowforc/low4/Ro0Fr1/n256/';
%namedir = '~/projects/pv/data_analysis/lowforc/low4/Ro1Fr0/n256/';
%namedir = '~/projects/pv/data_analysis/lowforc/low4/qg/qg256/bous500/';
%namedir = '~/projects/pv/data_analysis/lowforc/low4/qg/qg256/bous1000/';
%namedir = '~/projects/pv/data_analysis/lowforc/low4/qg/qg256/bous2000/';
%namedir = '~/projects/pv/data_analysis/lowforc/low4/qg/qg256/fcor2000_bous1000/';


%namedir = '~/projects/INCITE_runs/Intrepid/qg/';
%ame = 'n640_bous3000_all';
%kf=4;

namedir = '~/projects/INCITE_runs/Intrepid/Ro1Fr0/';
name = 'n640_fcor14bous3000_all';
kf=4;

%namedir = '~/projects/INCITE_runs/Intrepid/Ro0Fr1/';
%name = 'n640_fcor3000bous14_all';
%kf=4;

%namedir = '~/projects/INCITE_runs/Intrepid/bous_NSvisc/Ro1Fr0.01/';
%name = 'n512_Ro1Fr0.01_all';

%namedir = '~/projects/INCITE_runs/Intrepid/bous_NSvisc/Ro1Fr0.002/n1024_nu.2e-4/';
%name = 'n1024_Ro1Fr0.002_all';
%kf=4;

%namedir = '~/projects/INCITE_runs/Intrepid/bous_NSvisc/Ro1Fr0.002/n1024_nu.5e-5/';
%name = 'n1024_Ro1Fr0.002_all';

namedir = '~/projects/INCITE_runs/Intrepid/bous_NSvisc/Ro1Fr0.001/n1024_nu.5e-5/';
name = 'n1024_Ro1Fr0.001_all';
kf=4;

%namedir = '~/projects/INCITE_runs/Intrepid/bous_NSvisc/Ro0.002Fr1/n1024_nu.7e-5/';
%name = 'n1024_Ro0.002Fr1_all';
%kf=4;

%namedir = '~/projects/INCITE_runs/Intrepid/lowaspect_bous/LOWRES/sto_high_t4/';
%name = 'n512_d0.25_Ro0.05_all';
%kf=4;

%namedir = '~/projects/INCITE_runs/Intrepid/lowaspect_bous/n1024_d0.25_Ro0.05_nodamp/';
%name = 'n1024_d0.25_Ro0.05_all';
%kf=4;

%namedir = '~/projects/INCITE_runs/Intrepid/lowaspect_bous/n2048_d0.25_Ro0.05_nodamp/';
%name = 'n2048_d0.25_Ro0.05_0000.0000';
%kf=4;

asciprint = 0 % if == 1 print out the data to asci files

fid=fopen([namedir,name,'.spec2d'],'r','b');  %use 'b' for 
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
  spec2d_u=zeros([numkz,numkh]);  
  spec2d_v=zeros([numkz,numkh]);
  spec2d_w=zeros([numkz,numkh]);
  spec2d_t = zeros([numkz,numkh]);
  
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
        spec2d_u(kz,:) = s';
    elseif (ivar==2)
        spec2d_v(kz,:) = s';
    elseif (ivar==3)
        spec2d_w(kz,:) = s';
    elseif (ivar ==4)
        spec2d_t(kz,:) = s';
    end
    end
  end

  % kz-averaged 2d spectra (remaining array is function of kh)
  specuh = sum(spec2d_u,1);
  specvh = sum(spec2d_v,1);
  specwh = sum(spec2d_w,1);
  specth = sum(spec2d_t,1);

  %kh-averaged 1d spectra (remaining array is function of kz)
  specuz = sum(spec2d_u,2);
  specvz = sum(spec2d_v,2);
  specwz = sum(spec2d_w,2);
  spectz = sum(spec2d_t,2);

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
  loglog(kh,(spec2d_u(kzvals(i),:)'+spec2d_v(kzvals(i),:)').*(kh'.^expo)); hold on;%pause
  ylabel('E_h(k_h,k_z)')
  subplot(2,1,2); 
  loglog(kh,spec2d_t(kzvals(i),:)'.*(kh'.^expo));hold on;%pause
  ylabel('P(k_h,k_z)')
  figure(1); xlabel('k_h')
  end  
 

  %plot E_h(kh,kz) and P(kh,kz) as a function of kz for various kh
  figure(2);
  set(gca,'fontsize',16);
  subplot(2,1,1); hold off;
  subplot(2,1,2); hold off;
  for i = 1:length(khvals)
  subplot(2,1,1)
  loglog(kz,(spec2d_u(:,khvals(i))+spec2d_v(:,khvals(i))).*(kz'.^(expo)));hold on;%pause
  ylabel('E_h(k_h,k_z)')
  subplot(2,1,2)
  loglog(kz,spec2d_t(:,khvals(i)).*(kz'.^(expo)));hold on;%pause
  ylabel('P(k_h,k_z)')
  figure(2); xlabel('k_z')
  end
  

end

if (0)
  %plot the kz-averaged spectra
  figure(3)
  loglog53(numkh,specuh'+specvh','E_h(kh)',2.0,6);%hold on;
%  figure(4)
%  loglog53(numkh,specwh','E_z(kh)',2.0,6);%hold on;
  figure(5)
  loglog53(numkh,specth','P(kh)',2.0,6);%hold on;

  %plot the kh-averaged spectra
  figure(6)
  loglog53(numkz,specuz+specvz,'E_h(kz)',2.0,6);%hold on;
  figure(7)
  loglog53(numkz,specwz,'E_z(kz)',2.0,6);%hold on;
  figure(8)
  loglog53(numkz,spectz,'P(kz)',2.0,6);%hold on;
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
 
 %time average of spectra
 if(1)
 if (time > 1 & time < 5)
    if (j==0) 
        spec2d_t_ave = spec2d_t;
        spec2d_Eh_ave = spec2d_u + spec2d_v;
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
loglog(kz,spec2d_Eh_ave(:,khvals(i)).*((kz'.^expo)));hold on;pause;
ylabel('E_h(k_h,k_z)')
subplot(2,1,2)
axis([1 640/3 1e-10 1])
set(gca,'fontsize',17);
loglog(kz,spec2d_t_ave(:,khvals(i)).*((kz'.^expo)));hold on;pause;
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


if (forpaper == 1) 
 
figure(30);hold off;
for i = 1:length(khvals)     
axis([1 640/3 1e-10 10^5])
set(gca,'fontsize',17);
loglog(kz,10e4*spec2d_Eh_ave(:,khvals(i)).*((kz'.^expo)));hold on;pause;
loglog(kz,spec2d_t_ave(:,khvals(i)).*((kz'.^expo)),'r--');hold on;pause;
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
loglog(kh,specEhkh_ave.*((kh.^expo)));hold on;%pause
loglog(kh,specPkh_ave.*((kh.^expo)),'r--');hold on;%pause
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
loglog(kz,specEhkz_ave.*((kz'.^expo)));hold on;%pause
loglog(kz,specPkz_ave.*((kz'.^expo)),'r--');hold on;%pause
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