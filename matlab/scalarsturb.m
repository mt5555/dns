%
%########################################################################
%#  plot of DNS structure function output file
%########################################################################
%
%
clear all;
apply_fix=0;

range=0:50;
%range = 1
%range=0:.5:2;

fid=fopen('iso12_256_200.scalars-turb','r','b'); 


while (1) 
  [ns_e,count] = fread(fid,1,'float64');
  if (count~=1) break; end;

  time = fread(fid,1,'float64');
  data1 = fread(fid,[ns_e,1],'float64');
  data1=[time;data1];
  if (nscalars_e==0) 
    ints_e= data1;
  else
    ints_e= [ints_e,data1];
  end
  nscalars_e=nscalars_e+1;
  [nscalars_e,ns_e]
end
fclose(fid);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  the "expensive" scalars 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time_e=ints_e(1,:);
ux2=zeros([3,nscalars_e]);
ux3=ux2;
ux4=ux2;
for i=1:3
   ux2(i,:)=ints_e(i+1,:);    % < u1,1 ^ 2 >
   ux3(i,:)=ints_e(i+4,:);    % < u1,1 ^ 3 >
   ux4(i,:)=ints_e(i+7,:);    % < u1,1 ^ 4 >
end
Sww=ints_e(11,:);
if (ns_e>=13) 
  for i=1:3
    u2(i,:)=ints_e(i+11,:);              % < u^2 >
  end
end
if (ns_e>=14) 
    S2=ints_e(14,:);                     % <S2>
end


figure(6)
clf
hold on
%plot(time_e,Sww)

%plot(time_e,ux2(1,:),'k')
%plot(time_e,ux2(2,:),'b')
%plot(time_e,ux2(3,:),'g')

plot(time_e,ux3(1,:)./ux2(1,:).^(3/2),'k')
plot(time_e,ux3(2,:)./ux2(2,:).^(3/2),'b')
plot(time_e,ux3(3,:)./ux2(3,:).^(3/2),'g')

plot(time_e,ux4(1,:)./ux2(1,:).^2,'r')
plot(time_e,ux4(2,:)./ux2(2,:).^2,'b')
plot(time_e,ux4(3,:)./ux2(3,:).^2,'g')


hold off














