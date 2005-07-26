
name = 'temp0000.0000'
namedir = '../src/';

fid=endianopen([namedir,name,'.spec2d'],'r');
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

  disp(sprintf('time=%f  kz=%f  kh=%f',time,numkz,numkh));

  for kz=1:numkz
    spec2d(kz,:) = fread(fid,numkh,'float64') ;
  end
  
  loglog(spec2d(1,:))
  pause
end

