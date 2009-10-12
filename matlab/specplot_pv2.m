clear

name = 'n640_bous3000_all'
namedir = '~/INCITE_runs/Intrepid/qg/';


asciprint = 0 % if == 1 print out the data to asci files

fid=fopen([namedir,name,'.normpvspec'],'r','b');  %use 'b' for 
if (fid<0) 
  'Error opening file',[namedir,name]
  return
end

%read 0.5*|kz*theta|^2, 0.5*|kh*uh|^2 and potential enstrophy

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


if(1)

  %plot pv2(;,kh)/kztheta(:,kh) for various kh and pv2(kz,;)/khuh(kz,:) for various kz
  figure(1);hold off;
  subplot(2,1,1);
  kz = (1:numkz)-1
  kh = (1:numkh)-1
  xlabel('k_z')
  loglog(kz,sqrt(pv2_2d(:,3)./kztheta_2d(:,3)));hold on;pause
  loglog(kz,sqrt(pv2_2d(:,5)./kztheta_2d(:,5)));hold on;pause
  loglog(kz,sqrt(pv2_2d(:,11)./kztheta_2d(:,11)));hold on;pause
  loglog(kz,sqrt(pv2_2d(:,21)./kztheta_2d(:,21)));hold on;pause
  loglog(kz,sqrt(pv2_2d(:,31)./kztheta_2d(:,31)));hold on;pause
  loglog(kz,sqrt(pv2_2d(:,41)./kztheta_2d(:,41)));hold on;pause
  loglog(kz,sqrt(pv2_2d(:,51)./kztheta_2d(:,51)));hold on;pause
  loglog(kz,sqrt(pv2_2d(:,61)./kztheta_2d(:,61)));hold on;
  xlabel('k_z')
  ylabel('(q^2(k_z,k_h)/(k_z*\theta(k_z,k_h))^2)^{1/2}')
  

  subplot(2,1,2)
  xlabel('k_h') 
  loglog(kh,sqrt(pv2_2d(3,:)./kztheta_2d(3,:)));hold on;pause
  loglog(kh,sqrt(pv2_2d(5,:)./kztheta_2d(5,:)));hold on;pause
  loglog(kh,sqrt(pv2_2d(11,:)./kztheta_2d(11,:)));hold on;pause
  loglog(kh,sqrt(pv2_2d(21,:)./kztheta_2d(21,:)));hold on;pause
  loglog(kh,sqrt(pv2_2d(31,:)./kztheta_2d(31,:)));hold on;pause
  loglog(kh,sqrt(pv2_2d(41,:)./kztheta_2d(41,:)));hold on;pause
  loglog(kh,sqrt(pv2_2d(51,:)./kztheta_2d(51,:)));hold on;pause
  loglog(kh,sqrt(pv2_2d(61,:)./kztheta_2d(61,:)));hold on;
    xlabel('k_h') 
  
  %plot norm2spec(kz,:) for various kh/kz and kz/kh 
  figure(2);hold off;
  subplot(2,1,1)
  loglog(kz,sqrt(pv2_2d(:,3)./khuh_2d(:,3)));hold on;pause
  loglog(kz,sqrt(pv2_2d(:,5)./khuh_2d(:,5)));hold on;pause
  loglog(kz,sqrt(pv2_2d(:,11)./khuh_2d(:,11)));hold on;pause
  loglog(kz,sqrt(pv2_2d(:,21)./khuh_2d(:,21)));hold on;pause
  loglog(kz,sqrt(pv2_2d(:,31)./khuh_2d(:,31)));hold on;pause
  loglog(kz,sqrt(pv2_2d(:,41)./khuh_2d(:,41)));hold on;pause
  loglog(kz,sqrt(pv2_2d(:,51)./khuh_2d(:,51)));hold on;pause
  loglog(kz,sqrt(pv2_2d(:,61)./khuh_2d(:,61)));hold on;
  xlabel('k_z');
  
  subplot(2,1,2)
  loglog(kh,sqrt(pv2_2d(3,:)./khuh_2d(3,:)));hold on;pause
  loglog(kh,sqrt(pv2_2d(5,:)./khuh_2d(5,:)));hold on;pause
  loglog(kh,sqrt(pv2_2d(11,:)./khuh_2d(11,:)));hold on;pause
  loglog(kh,sqrt(pv2_2d(21,:)./khuh_2d(21,:)));hold on;pause
  loglog(kh,sqrt(pv2_2d(31,:)./khuh_2d(31,:)));hold on;pause
  loglog(kh,sqrt(pv2_2d(41,:)./khuh_2d(41,:)));hold on;pause
  loglog(kh,sqrt(pv2_2d(51,:)./khuh_2d(51,:)));hold on;pause
  loglog(kh,sqrt(pv2_2d(61,:)./khuh_2d(61,:)));hold on;
  xlabel('k_h');
  
  
  
  
end



