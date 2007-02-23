%
%read the spectra of potential enstrophy from .pv2spec files. That
%file contains ALL potential enstrophy spectra, shell-averaged (1 dim),
%plane-averaged (1 dim) and annulus averaged (2 dim).
%

clear
%name = 'all'
%name = 'qg256hyper_all32'
%name = 'qg64hyper_all80'
name = 'qg256hyper0000.0000'
%epsilon=.41;
%CK=1.5*epsilon^(2/3);
%namedir = '/home/wingate/data1/Rotation/r16c/';
namedir = '~/projects/pv/data_analysis/lowforc/low4/qg256/bous1000/hyper_nu15/';
%namedir = '~/projects/pv/data_analysis/lowforc/low4/qg64/sto_high_4/hyper_nu/bous500/';



fid=fopen([namedir,name,'.pv2spec'],'r');
if (fid<0) 
  'Error opening file',[namedir,name]
  return
end


time=0;
while (1)
  [time,count]=fread(fid,1,'float64');


  if (count ~= 1) 
    disp('EOF reached.  stopping')
    break 
  end
  if (time<0 | time>1000) break; end;
  numkr=fread(fid,1,'float64');
  numkx=fread(fid,1,'float64');
  numky=fread(fid,1,'float64');
  numkz=fread(fid,1,'float64');
  numkh=fread(fid,1,'float64');
  q2spec_r=zeros([numkr]);  
  q2spec_x=zeros([numkx]);
  q2spec_y=zeros([numky]);
  q2spec_z=zeros([numkz]);
  q2spec_2d=zeros([numkz,numkh]);
  
  disp(sprintf('time=%f  kz=%f  kh=%f',time,numkz,numkh));


  q2spec_r = fread(fid,numkr,'float64') ;
  q2spec_x = fread(fid,numkx,'float64') ;
  q2spec_y = fread(fid,numky,'float64') ;
  q2spec_z = fread(fid,numkz,'float64') ; 
  
  
  for kz=1:numkz
    [q2spec_2d(kz,:),count] = fread(fid,numkh,'float64') ;
    if (count~=numkh)
      disp('Error: error reading file')
    end
  end

  % kz-averaged 2d spectrum (becomes 1-dimensional)
  q2spec_2d_ave = sum(q2spec_2d,1)/numkz;

 
  %plot the shell-averaged and plane-averaged spectra
  figure(1)
  loglog53(numkr,q2spec_r','',2.0,6);hold on;
  pause
  loglog53(numkx,q2spec_x','',2.0,6);hold on;
  pause
  loglog53(numky,q2spec_y','',2.0,6);hold on;
  pause
  loglog53(numkz,q2spec_z','',2.0,6);hold on;
  
  %plot the kz_averaged annulus spectrum
  figure(2)
  loglog53(numkh,q2spec_2d_ave','',2.0,6);hold on;
  
  pause
end


