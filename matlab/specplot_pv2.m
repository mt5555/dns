clear

name = 'n640_bous3000_0006.0000'
namedir = '~/INCITE_runs/Intrepid/qg/';


asciprint = 0 % if == 1 print out the data to asci files

fid=fopen([namedir,name,'.normpvspec'],'r','b');  %use 'b' for 
if (fid<0) 
  'Error opening file',[namedir,name]
  return
end

  %read norm1spec = 0.5*|kz*theta|^2, and norm2psec = 0.5*|kh*uh|^2

time=0;
j = 0; %(count for time average)
while (time < 999)
  [time,count]=fread(fid,1,'float64');
  if (count ~= 1) 
    disp('EOF reached.  stopping')
    break 
  end
  if (time<0 | time>1000) break; end;
  numkz=fread(fid,1,'float64');
  numkh=fread(fid,1,'float64');
    
kztheta_2d = zeros([numkz,numkh]);

  disp(sprintf('time=%f  kz=%f  kh=%f',time,numkz,numkh));
for kz=1:numkz        
  [s,count] = fread(fid,numkh,'float64') ;
if (count~=numkh)
     disp('Error: error reading file')
     count
size(s)
     end
        kztheta_2d(kz,:) = s';
end



khuh_2d = zeros([numkz,numkh]);

for kz=1:numkz        
  [s,count] = fread(fid,numkh,'float64') ;
if (count~=numkh)
     disp('Error: error reading file')
     count
size(s)
     end
        khuh_2d(kz,:) = s';
end


pv2_2d = zeros([numkz,numkh]);

for kz=1:numkz        
  [s,count] = fread(fid,numkh,'float64') ;
if (count~=numkh)
     disp('Error: error reading file')
     count
size(s)
     end
        pv2_2d(kz,:) = s';
end

end

  % kz-averaged 2d spectra (remaining array is function of kh)
%  n1spec = sum(norm1spec_r_2d,1);
  
  %kh-averaged 1d spectra (remaining array is function of kz)
%  n2spec = sum(norm2spec_r_2d,2);

if(1)

  %plot norm1spec(:,kh) for various kh and norm1spec(kz,:) for various kh
  figure(1);hold off;
  subplot(2,1,1);
  kz = (1:numkz)-1
  kh = (1:numkh)-1
  xlabel('k_z')
  loglog(kz,pv2_2d(:,3)./kztheta_2d(:,3));hold on;pause
  loglog(kz,pv2_2d(:,5)./kztheta_2d(:,5));hold on;pause
  loglog(kz,pv2_2d(:,11)./kztheta_2d(:,11));hold on;pause
  loglog(kz,pv2_2d(:,21)./kztheta_2d(:,21));hold on;pause
  loglog(kz,pv2_2d(:,31)./kztheta_2d(:,31));hold on;pause
  loglog(kz,pv2_2d(:,41)./kztheta_2d(:,41));hold on;pause
  loglog(kz,pv2_2d(:,51)./kztheta_2d(:,51));hold on;pause
  loglog(kz,pv2_2d(:,61)./kztheta_2d(:,61));hold on;
  xlabel('k_z')
  
  subplot(2,1,2)
  xlabel('k_h') 
  loglog(kh,pv2_2d(3,:)./kztheta_2d(3,:));hold on;pause
  loglog(kh,pv2_2d(5,:)./kztheta_2d(5,:));hold on;pause
  loglog(kh,pv2_2d(11,:)./kztheta_2d(11,:));hold on;pause
  loglog(kh,pv2_2d(21,:)./kztheta_2d(21,:));hold on;pause
  loglog(kh,pv2_2d(31,:)./kztheta_2d(31,:));hold on;pause
  loglog(kh,pv2_2d(41,:)./kztheta_2d(41,:));hold on;pause
  loglog(kh,pv2_2d(51,:)./kztheta_2d(51,:));hold on;pause
  loglog(kh,pv2_2d(61,:)./kztheta_2d(61,:));hold on;
    xlabel('k_h') 
  
  %plot norm2spec(kz,:) for various kh/kz and kz/kh 
  figure(2);hold off;
  subplot(2,1,1)
  loglog(kz./2,norm2_2d(:,3));hold on;pause
  loglog(kz./4,norm2_2d(:,5));hold on;pause
  loglog(kz./10,norm2_2d(:,11));hold on;pause
  loglog(kz./20,norm2_2d(:,21));hold on;pause
  loglog(kz./30,norm2_2d(:,31));hold on;pause
  loglog(kz./40,norm2_2d(:,41));hold on;pause
  loglog(kz./50,norm2_2d(:,51));hold on;pause
  loglog(kz./60,norm2_2d(:,61));hold on;
  xlabel('\kappa_z/k_h');
  
  subplot(2,1,2)
  loglog(2./kh,norm2_2d(3,:)');hold on;pause
  loglog(4./kh,norm2_2d(5,:)');hold on;pause
  loglog(10./kh,norm2_2d(11,:)');hold on;pause
  loglog(20./kh,norm2_2d(21,:)');hold on;pause
  loglog(30./kh,norm2_2d(31,:)');hold on;pause
  loglog(40./kh,norm2_2d(41,:)');hold on;pause
  loglog(50./kh,norm2_2d(51,:)');hold on;pause
  loglog(60./kh,norm2_2d(61,:)');hold on;
  xlabel('\kappa_z/k_h');
  
  
  
  
end



