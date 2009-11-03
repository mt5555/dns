clear

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


%name = 'n640_fcor3000bous14_all';
%name = 'n640_fcor14bous3000_all';
%name = 'n640_bous3000_all';
name = 'n512_Ro1Fr0.01_all';

%namedir = '~/projects/INCITE_runs/Intrepid/qg/';
%namedir = '~/projects/INCITE_runs/Intrepid/Ro1Fr0/';
%namedir = '~/projects/INCITE_runs/Intrepid/Ro0Fr1/';
namedir = '~/projects/INCITE_runs/Intrepid/bous_NSvisc/';

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

if(1)
  %plot E_h(0,kh)
  kz = (1:numkz)-1;
  kh = (1:numkh)-1;
  figure(1);hold off;
  expo = 4.5 ; %0 for no compensation
  kzvals = [1,3,5,11,21,31];
  for i = 1:length(kzvals)
  loglog(kh,(spec2d_u(kzvals(i),:)'+spec2d_v(kzvals(i),:)').*(kh'.^(expo)));hold on;%pause
  end  
 

  %plot P(kz,0)
  figure(2);hold off;
  expo = 5;
  khvals = [1,3,5,11,21,31];
  for i = 1:length(khvals)
  loglog(kz,spec2d_t(:,khvals(i)).*(kz'.^(expo)));hold on;%pause
  end
  ylabel('P(:,[1,3,5,11,21,31,41,51])');  
  xlabel('kz');

end

if (1)
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
%  figure(7)
%  loglog53(numkz,specwz,'E_z(kz)',2.0,6);%hold on;
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
 if (time > 3.0 & time < 6.5)
    if (j==0) 
        spec2d_t_ave = spec2d_t;
        spec2d_Eh_ave = spec2d_u + spec2d_v;
        spectz_ave = spectz;
        j = j+1;
    else
        spec2d_t_ave = spec2d_t_ave + spec2d_t;
        spec2d_Eh_ave = spec2d_Eh_ave + spec2d_u + spec2d_v;
        spectz_ave = spectz_ave + spectz;
        
        j = j + 1;
    end
 end
 end 
end
spectz_ave = spectz_ave/j;
spec2d_t_ave = spec2d_t_ave/j;
spec2d_Eh_ave = spec2d_Eh_ave/j;

%plot time avg spectra

figure(10);hold off;
expo = 5;
for i = 1:length(khvals) 
loglog(kz,spec2d_t_ave(:,khvals(i)).*(kz'.^expo));hold on;%pause
end 
title('Time avg P(kh,kz)')
   
figure(11);hold off;
expo = 4.5;
for i = 1:length(kzvals) 
loglog(kh,spec2d_Eh_ave(kzvals(i),:).*(kh.^expo));hold on;%pause
end 
title('Time avg E_h(kh,kz)')

%normalize x-axis by ratio kz or kh as necessary
figure(12)
for kh = 10:10:50
loglog(((1:numkz)-1)/kh,spec2d_t_ave(:,kh+1),'',2.0,6);hold on;%pause
xlabel('$kz/kh$')
ylabel('$P(k_h,k_z)$')
title('P(kh,kz) vs. kz/kh')
end

figure(13)
for kz = 10:10:50
loglog(((1:numkh)-1)/kz,spec2d_Eh_ave(kz+1,:)','',2.0,6);hold on;%pause
xlabel('$kh/kz$')
ylabel('$E_h(k_h,k_z)$')
title('Eh(kh,kz) vs. kh/kz')

end


