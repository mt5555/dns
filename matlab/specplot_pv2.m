clear

name = 'n640_bous3000_all'
namedir = '~/projects/INCITE_runs/Intrepid/qg/';


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

%disp(sprintf('time=%f  k=%f  kx=%f  ky=%f  kz=%f  kh=%f',time,numk,numkx,numky,numkz,numkh));    
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

% normalized quantities
qnorm1 = sqrt(pv2_2d./kztheta_2d);
qnorm2 = sqrt(pv2_2d./khuh_2d);

% time accumulate   
   if(1)
 if (time > 3.0 & time < 6.5)
    if (j==0) 
        qnorm1_ave = sqrt(pv2_2d./kztheta_2d);
        qnorm2_ave = sqrt(pv2_2d./khuh_2d);
        j = j+1;
    else
        qnorm1_ave = qnorm1_ave + sqrt(pv2_2d./kztheta_2d);
        qnorm2_ave = qnorm2_ave + sqrt(pv2_2d./khuh_2d);
        j = j+1;
        
    end
 end
 end 
end

%time average
qnorm1_ave = qnorm1_ave/j;
qnorm2_ave = qnorm2_ave/j;

end

if(1)

  %plot pv2(;,kh)/kztheta(:,kh) for various kh and pv2(kz,;)/khuh(kz,:) for various kz
  figure(1);hold off;
  subplot(1,1,1);
  kz = (1:numkz)-1;
  kh = (1:numkh)-1;
  xlabel('k_z');
  khvals = [1,2,3,4,5,11,21]
  for i = 1:length(khvals)
  loglog(kz,qnorm1(:,khvals(i)),'-','Linewidth',1.5);hold on;pause;
  end
  set(gca,'fontsize',16);    
  xlabel('k_z')
  ylabel('(|q(k)|^2/|k_z \theta(k)|^2)^{1/2}')
  

  %plot q2/khuh(kz,:) for various kz 
   figure(2);hold off;
  subplot(1,1,1)
  kzvals = [1,2,3,4,5,11,21]
  for i = 1:length(kzvals)
  loglog(kh,qnorm2(kzvals(i),:),'-','Linewidth',1.5);hold on;pause;
  end
    set(gca,'fontsize',16);  
    xlabel('k_h');
  ylabel('(|q(k)|^2/|k_h u_h(k)|^2)^{1/2}')  
  
  % time average norm spectra
  figure(3);hold off
  for i = 1:length(khvals)
    loglog(kz,qnorm1_ave(:,khvals(i)),'-','Linewidth',1.5);hold on;pause;
  end
    set(gca,'fontsize',16);
  xlabel('k_z')
  ylabel('(|q(k)|^2/|k_z \theta(k)|^2)^{1/2}')
  
  figure(4);hold off
  for i = 1:length(kzvals)
    loglog(kh,qnorm2_ave(kzvals(i),:),'-','Linewidth',1.5);hold on;pause;
  end
  set(gca,'fontsize',16);  
  xlabel('k_h');
  ylabel('(|q(k)|^2/|k_h u_h(k)|^2)^{1/2}');  
%pv2spectra
   figure(5);hold off;
   k = (1:numk)-1;
   loglog(k,q2_r);  
  

end




