clear
name = 'r160000.0000'
epsilon=.41;
CK=1.5*epsilon^(2/3);
namedir = '/home/mataylo/';

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
  numkh=fread(fid,1,'float64');
  numkz=fread(fid,1,'float64');
  spec2d=zeros([numkz,numkh]);  

  disp(sprintf('time=%f  kz=%f  kh=%f',time,numkz,numkh));

  for kz=1:numkz
    [s,count] = fread(fid,numkh,'float64') ;
    if (count~=numkh)
      disp('Error: error reading file')
      count
      size(s)
    end
    spec2d(kz,:) = s';
  end
%  for kz=1:numkz
%    for kh = 1,numkh
%      spec2d(kz,kh) = fread(fid,numkh+1,'float64') ;
%    end
%  end

 
  figure(1)
  loglog53(numkh,spec2d(1,:)','',2.0,6)

  hold on;
  spec = sum(spec2d,1)/numkz;
  loglog53(numkh,spec','E0(kh) and E(kh)',2.0,6)

  
  pause
end

