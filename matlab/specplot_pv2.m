clear

name = 'n640_bous3000_all'
namedir = '~/INCITE_runs/Intrepid/qg/';


asciprint = 0 % if == 1 print out the data to asci files

fid=fopen([namedir,name,'.normspec'],'r','b');  %use 'b' for 
if (fid<0) 
  'Error opening file',[namedir,name]
  return
end

  %read norm1spec = q/kz*theta, and norm2psec = q/kh*uh

time=0;
j = 0; %(count for time average)
while (time < 999)
  [time,count]=fread(fid,1,'float64');
  if (count ~= 1) 
    disp('EOF reached.  stopping')
    break 
  end
  if (time<0 | time>1000) break; end;
%  numvar=fread(fid,1,'float64');
  numk=fread(fid,1,'float64');
  numkx=fread(fid,1,'float64');
  numky=fread(fid,1,'float64');
  numkz=fread(fid,1,'float64');
  numkh=fread(fid,1,'float64');
  
  norm1spec_r=zeros([numk]);  
  norm1spec_x=zeros([numkx]);
  norm1spec_y=zeros([numky);
  norm1spec_z = zeros([numkz]);
norm1_2d = zeros([numkz,numkh]);

norm1spec_r = fread(fid,numk,'float64');
norm1spec_x = fread(fid,numkx,'float64');
norm1spec_y = fread(fid,numky,'float64');
norm1spec_z = fread(fid,numkz,'float64');

  disp(sprintf('time=%f  kz=%f  kh=%f',time,numkz,numkh));
for kz=1:numkz        
  [s,count] = fread(fid,numkh,'float64') ;
if (count~=numkh)
     disp('Error: error reading file')
     count
size(s)
     end
        norm1_2d(kz,:) = s';
end

  
  norm2spec_r=zeros([numk]);  
  norm2spec_x=zeros([numkx]);
  norm2spec_y=zeros([numky]);
  norm2spec_z = zeros([numkz]);
  norm2_2d = zeros([numkz,numkh]);

norm2spec_r = fread(fid,numk,'float64');
norm2spec_x = fread(fid,numkx,'float64');
norm2spec_y = fread(fid,numky,'float64');
norm2spec_z = fread(fid,numkz,'float64');


for kz=1:numkz        
  [s,count] = fread(fid,numkh,'float64') ;
if (count~=numkh)
     disp('Error: error reading file')
     count
size(s)
     end
        norm2_2d(kz,:) = s';
end



  % kz-averaged 2d spectra (remaining array is function of kh)
%  n1spec = sum(norm1spec_r_2d,1);
  
  %kh-averaged 1d spectra (remaining array is function of kz)
%  n2spec = sum(norm2spec_r_2d,2);

if(1)

  %plot norm1spec(:,kh) for various kh 
  figure(1);hold off;
  loglog53(numkh,norm1_2d(:,1)','',2.0,6);hold on;%pause
  loglog53(numkh,norm1_2d(:,3)','',2.0,6);hold on;%pause
  loglog53(numkh,norm1_2d(:,5)','',2.0,6);hold on;%pause
  loglog53(numkh,norm1_2d(:,11)','',2.0,6);hold on;%pause
  loglog53(numkh,norm1_2d(:,21)','',2.0,6);hold on;%pause
  loglog53(numkh,norm1_2d(:,31)','',2.0,6);hold on;%pause
  loglog53(numkh,norm1_2d(:,41)','',2.0,6);hold on;%pause
  loglog53(numkh,norm1_2d(:,51)','',2.0,6);hold on;%pause
  loglog53(numkh,norm1_2d(:,61)','|q/kz*theta| as function of kz for various kh',2.0,6);hold on;


  %plot norm2spec(kz,:) for various kz 
  figure(2);hold off;
  loglog53(numkh,norm2_2d(1,:)','',2.0,6);hold on;%pause
  loglog53(numkh,norm2_2d(3,:)','',2.0,6);hold on;%pause
  loglog53(numkh,norm2_2d(5,:)','',2.0,6);hold on;%pause
  loglog53(numkh,norm2_2d(11,:)','',2.0,6);hold on;%pause
  loglog53(numkh,norm2_2d(21,:)','',2.0,6);hold on;%pause
  loglog53(numkh,norm2_2d(31,:)','',2.0,6);hold on;%pause
  loglog53(numkh,norm2_2d(41,:)','',2.0,6);hold on;%pause
  loglog53(numkh,norm2_2d(51,:)','',2.0,6);hold on;%pause
  loglog53(numkh,norm2_2d(61,:)','|q/kh*uh| as function of kh for various kz',2.0,6);hold on;

end



