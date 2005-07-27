
name = 'temp0000.0000'
epsilon=.41;
CK=1.5*epsilon^(2/3);
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

  loglog53(numkh,spec2d(1,:)','',2.0,6)

  hold on;
  spec = sum(spec2d,1);
  scale =  (spec2d(1,5)/spec(1,5))
  spec = spec*scale;
  loglog53(numkh,spec','nz*E0(kh) and E(kh)',2.0,6)

  
  pause
end

