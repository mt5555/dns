clear
%name = 'all'
name = 'qg512hyper_all';
%name = 'qg256hyper_all';
%name = 'qg64hyper_all';
%name = 'n256_f1000n5_all';
%name = 'n256_f2000n5_all';
%name = 'n256_f5n1000_all';
%name = 'n256_f5n2000_all';
%name = 'n256high_f2000n5_all';
%name = 'n256high_f1000n5_all'
%name = 'n256high_f5n2000_all';
%name = 'n256high_f5n1000_all'

epsilon=.41;
%CK=1.5*epsilon^(2/3);
%namedir = '/home/wingate/data1/Rotation/r16c/';
%namedir = '~/projects/pv/data_analysis/lowforc/low4/qg256/bous1000/hyper_nu2.5/';
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
namedir = '~/projects/pv/data_analysis/lowforc/low4/qg/qg512/bous2000/';
asciprint = 0 % if == 1 print out the data to asci files

fid=fopen([namedir,name,'.spec2d'],'r');
if (fid<0) 
  'Error opening file',[namedir,name]
  return
end


time=0;
j = 0; %(count for time average)
while (time < 4.5)
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
  specu = sum(spec2d_u,1);
  specv = sum(spec2d_v,1);
  specw = sum(spec2d_w,1);
  spect = sum(spec2d_t,1);

  %kh-averaged 2d spectra (remaining array is function of kz)
  specuz = sum(spec2d_u,2);
  specvz = sum(spec2d_v,2);
  specwz = sum(spec2d_w,2);
  spectz = sum(spec2d_t,2);

if(1)
  %plot E_h(0,kh)
  figure(1);hold off;
  loglog53(numkh,spec2d_u(1,:)'+spec2d_v(1,:)','',2.0,6);hold on;%pause
  loglog53(numkh,spec2d_u(10,:)'+spec2d_v(11,:)','',2.0,6);hold on;%pause
  loglog53(numkh,spec2d_u(10,:)'+spec2d_v(21,:)','',2.0,6);hold on;%pause
  loglog53(numkh,spec2d_u(10,:)'+spec2d_v(31,:)','',2.0,6);hold on;%pause
  loglog53(numkh,spec2d_u(50,:)'+spec2d_v(41,:)','E_h((0, 10, 20, 30, 40),kh)',2.0,6);hold on;

  %plot P(kz,0)
  figure(2);hold off;
  loglog53(numkz,spec2d_t(:,1),'',2.0,6);hold on;%pause
  loglog53(numkz,spec2d_t(:,11),'',2.0,6);hold on;%pause
  loglog53(numkz,spec2d_t(:,21),'',2.0,6);hold on;%pause
  loglog53(numkz,spec2d_t(:,31),'',2.0,6);hold on;%pause
  loglog53(numkz,spec2d_t(:,41),'P(kz,(0,10,20,30,40,50))',2.0,6);hold on;
   axis([1 128 1e-12 1e-4]) ;
end

if (0)
  %plot the kz-averaged spectra
  figure(3)
  loglog53(numkh,specu'+specv','E_h(kh)',2.0,6);%hold on;
%  figure(4)
%  loglog53(numkh,specw','E_z(kh)',2.0,6);%hold on;
%  figure(5)
%  loglog53(numkh,spect','P(kh)',2.0,6);%hold on;

  %plot the kh-averaged spectra
%  figure(6)
%  loglog53(numkz,specuz+specvz,'E_h(kz)',2.0,6);%hold on;
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
 disp('pause') 
 pause
 
 %time average of spectra
 if (time > 1 & time < 2.0)
    if (j==0) 
        spec2d_t_ave = spec2d_t;
        spectz_ave = spectz;
        j = j+1;
    else
        spec2d_t_ave = spec2d_t_ave + spec2d_t;
        spectz_ave = spectz_ave + spectz;
        j = j + 1;
    end
 end
end
spectz_ave = spectz_ave/j;
spec2d_t_ave = spec2d_t_ave/j;

%plot time avg spectra
if(0)
figure(9)
loglog53(numkz,spectz_ave,'time-average P(kz)',2.0,6);

figure(10);hold off;
  loglog53(numkz,spec2d_t_ave(:,1),'',2.0,6);hold on;%pause
  loglog53(numkz,spec2d_t_ave(:,11),'',2.0,6);hold on;%pause
  loglog53(numkz,spec2d_t_ave(:,21),'',2.0,6);hold on;%pause
  loglog53(numkz,spec2d_t_ave(:,31),'',2.0,6);hold on;%pause
  loglog53(numkz,spec2d_t_ave(:,41),'P(kz,(0,10,20,30,40,50))',2.0,6);hold on;
   axis([1 128 1e-12 1e-4]) ;
end

