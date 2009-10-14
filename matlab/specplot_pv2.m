clear

name = 'n640_bous3000_all'
namedir = '~/INCITE_runs/Intrepid/qg/';


asciprint = 0 % if == 1 print out the data to asci files

fid0=fopen([namedir,name,'.pv2spec'],'r','b');  %use 'b' for 
if (fid0<0) 
  'Error opening file',[namedir,name]
end

fid=fopen([namedir,name,'.normpvspec'],'r','b');  %use 'b' for 
if (fid<0) 
  'Error opening file',[namedir,name]
  return
end

    if (fid0>0)
time=0;
j = 0; %(count for time average)
while (time < 999)
  [time,count]=fread(fid0,1,'float64');
  if (count ~= 1) 
    disp('EOF reached.  stopping')
    break 
  end
  if (time<0 | time>1000) break; end;
numk= fread(fid0,1,'float64');
  numkx=fread(fid0,1,'float64');
  numky=fread(fid0,1,'float64');
  numkz=fread(fid0,1,'float64');
  numkh=fread(fid0,1,'float64');

    
  q2_r = fread(fid0,numk,'float64') ;
q2_x = fread(fid0,numkx,'float64');
q2_y = fread(fid0,numky,'float64');
q2_z = fread(fid0,numkz,'float64');

end

end
%read 0.5*|kz*theta|^2, 0.5*|kh*uh|^2 and potential enstrophy

  if (fid>0)
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

end

if(1)

  %plot pv2(;,kh)/kztheta(:,kh) for various kh and pv2(kz,;)/khuh(kz,:) for various kz
  figure(1);hold off;
  subplot(2,1,1);
  kz = (1:numkz)-1;
  kh = (1:numkh)-1;
  xlabel('k_z')
  axis([1 1000 2000 5000]);
  semilogx(kz,sqrt(pv2_2d(:,3)./kztheta_2d(:,3)));hold on;pause
  semilogx(kz,sqrt(pv2_2d(:,5)./kztheta_2d(:,5)));hold on;pause
  semilogx(kz,sqrt(pv2_2d(:,11)./kztheta_2d(:,11)));hold on;pause
  semilogx(kz,sqrt(pv2_2d(:,21)./kztheta_2d(:,21)));hold on;pause
  semilogx(kz,sqrt(pv2_2d(:,31)./kztheta_2d(:,31)));hold on;pause
  semilogx(kz,sqrt(pv2_2d(:,41)./kztheta_2d(:,41)));hold on;pause
  semilogx(kz,sqrt(pv2_2d(:,51)./kztheta_2d(:,51)));hold on;pause
  semilogx(kz,sqrt(pv2_2d(:,61)./kztheta_2d(:,61)));hold on;
  xlabel('k_z')
  ylabel('(q^2(k_z,k_h)/(k_z \theta(k_z,k_h))^2)^{1/2}')
  

  subplot(2,1,2)
  xlabel('k_h') 
  axis([1 1000 2000 5000]);
  semilogx(kh,sqrt(pv2_2d(3,:)./kztheta_2d(3,:)));hold on;pause
  semilogx(kh,sqrt(pv2_2d(5,:)./kztheta_2d(5,:)));hold on;pause
  semilogx(kh,sqrt(pv2_2d(11,:)./kztheta_2d(11,:)));hold on;pause
  semilogx(kh,sqrt(pv2_2d(21,:)./kztheta_2d(21,:)));hold on;pause
  semilogx(kh,sqrt(pv2_2d(31,:)./kztheta_2d(31,:)));hold on;pause
  semilogx(kh,sqrt(pv2_2d(41,:)./kztheta_2d(41,:)));hold on;pause
  semilogx(kh,sqrt(pv2_2d(51,:)./kztheta_2d(51,:)));hold on;pause
  semilogx(kh,sqrt(pv2_2d(61,:)./kztheta_2d(61,:)));hold on;
    xlabel('k_h') 
    ylabel('(q^2(k_z,k_h)/(k_z \theta(k_z,k_h))^2)^{1/2}')

  %plot norm2spec(kz,:) for various kh/kz and kz/kh 
  figure(2);hold off;
  subplot(2,1,1)
  semilogx(kz,sqrt(pv2_2d(:,3)./khuh_2d(:,3)));hold on;pause
  semilogx(kz,sqrt(pv2_2d(:,5)./khuh_2d(:,5)));hold on;pause
  semilogx(kz,sqrt(pv2_2d(:,11)./khuh_2d(:,11)));hold on;pause
  semilogx(kz,sqrt(pv2_2d(:,21)./khuh_2d(:,21)));hold on;pause
  semilogx(kz,sqrt(pv2_2d(:,31)./khuh_2d(:,31)));hold on;pause
  semilogx(kz,sqrt(pv2_2d(:,41)./khuh_2d(:,41)));hold on;pause
  semilogx(kz,sqrt(pv2_2d(:,51)./khuh_2d(:,51)));hold on;pause
  semilogx(kz,sqrt(pv2_2d(:,61)./khuh_2d(:,61)));hold on;
  xlabel('k_z');
  ylabel('(q^2(k_z,k_h)/(k_h u_h(k_z,k_h))^2)^{1/2}')
  
  subplot(2,1,2)
  semilogx(kh,sqrt(pv2_2d(3,:)./khuh_2d(3,:)));hold on;pause
  semilogx(kh,sqrt(pv2_2d(5,:)./khuh_2d(5,:)));hold on;pause
  semilogx(kh,sqrt(pv2_2d(11,:)./khuh_2d(11,:)));hold on;pause
  semilogx(kh,sqrt(pv2_2d(21,:)./khuh_2d(21,:)));hold on;pause
  semilogx(kh,sqrt(pv2_2d(31,:)./khuh_2d(31,:)));hold on;pause
  semilogx(kh,sqrt(pv2_2d(41,:)./khuh_2d(41,:)));hold on;pause
  semilogx(kh,sqrt(pv2_2d(51,:)./khuh_2d(51,:)));hold on;pause
  semilogx(kh,sqrt(pv2_2d(61,:)./khuh_2d(61,:)));hold on;
  xlabel('k_h');
  ylabel('(q^2(k_z,k_h)/(k_h u_h(k_z,k_h))^2)^{1/2}')  
  
%pv2spectra
   figure(3);hold off;
   k = (1:numk)-1;
   loglog(k,q2_r);  
  
end



