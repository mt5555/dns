%
%########################################################################
%#  plot of DNS structure function output file
%########################################################################
%
%
clear all;

range=0:50;
%range = 1
%range=0:.5:2;

fid=fopen('test64.scalars','r');

nscalars=0;
while (1) 
  [ni,count]=fread(fid,1,'float64');
  if (count~=1) break; end;
  nints=ni;
  ns=fread(fid,1,'float64');
  
  data1=fread(fid,[nints,ns],'float64');
  data2=fread(fid,[nints,ns],'float64');
  if (nscalars==0) 
    ints=data1;
    maxs=data2;
  else
    ints=[ints,data1];
    maxs=[maxs,data2];
  end
  nscalars=nscalars+ns;

  
end

disp(sprintf('nints=%i  total scalars read=%i',nints,nscalars))
fclose(fid);


ke=ints(1,:);
ke_diss_f=ints(2,:);
ke_diss_d=ints(3,:);
vor=ints(4,:);
hel=ints(5,:);
ke_diss_tot=ints(6,:);

maxU=maxs(1,:);
maxV=maxs(2,:);
maxW=maxs(3,:);
maxvor=maxs(5,:);
timeU=maxs(6,:);
timeDU=maxs(7,:);


disp(sprintf('max vor = %e',max(vor)));

figure(5)
clf
hold on
plot(timeU,ke)
plot(timeDU(2:nscalars),ke_diss_tot(2:nscalars),'r')
plot(timeDU,ke_diss_f+ke_diss_d,'k')
%plot(timeDU,ke_diss_f,'k')
%plot(timeDU,ke_diss_d,'k')
plot(timeDU,hel,'g')
title('KE: blue,    d(KE)/dt: black & red,    hel: green');
hold off


















