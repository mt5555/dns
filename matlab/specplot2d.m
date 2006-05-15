clear
name = 'all'
epsilon=.41;
CK=1.5*epsilon^(2/3);
namedir = '/home/wingate/data1/Rotation/r16c/';

asciprint = 1 % if == 1 print out the data to asci files

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

 
  figure(1)
  loglog53(numkh,spec2d(1,:)','',2.0,6)

  hold on;
  spec = sum(spec2d,1)/numkz;
  loglog53(numkh,spec','E0(kh) and E(kh)',2.0,6)

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


