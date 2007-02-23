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

asciprint = 0 % if == 1 print out the data to asci files

fid=fopen([namedir,name,'.spec2d'],'r');
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

  % kz-averaged 2d spectra
  specu = sum(spec2d_u,1)/numkz;
  specv = sum(spec2d_v,1)/numkz;
  specw = sum(spec2d_w,1)/numkz;
  spect = sum(spec2d_t,1)/numkz;

  %plot the kz-averaged spectra
  figure(1)
  loglog53(numkh,specu'+specv','E_h(kh)',2.0,6);hold on;
  figure(2)
  loglog53(numkh,specw','E_z(kh)',2.0,6);hold on;
  figure(3)
  loglog53(numkh,spect','E_\theta(kh)',2.0,6);hold on;
  
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
  
  pause
end


