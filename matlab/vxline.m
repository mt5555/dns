%
%########################################################################
%#  plotting the .vxline data
%########################################################################
clear;

%ts=input('time=? ');

%name='/home/taylorm/ccs/dns/src/vxpair/vx4096d';
name='/home/taylorm/ccs/dns/src/vxpair/vx4500a';
times=[0:.1:180];

name='/scratch2/taylorm/vx12288b/vx12288b';
times=[0:.1:180];



ccol=[ 'b','g','r','c','m','y', 'b','g','r','c','m','y' ];  
i=find(name=='/');
i=i(length(i))
pname=name(i+1:length(name));

timev=[];
wlinev=zeros([10,1]);
k=0;


for t=times
  tstr=sprintf('%10.4f',t+10000);
  fname=[name,tstr(2:10),'.vxline'];
  fid=endianopen(fname,'r');

  if (fid>=0) 
     disp(fname)

     [nline,count]=fread(fid,1,'float64');
     if (count~=1) break; end;
     time=fread(fid,1,'float64');
     xc=fread(fid,1,'float64');
     yc=fread(fid,1,'float64');
     yline=fread(fid,nline,'float64');
     wline=fread(fid,nline,'float64');
     fclose(fid);
     k=k+1;
     timev(k)=time;
     wlinev(1:9,k)=wline(1:9) ;
     xcv(k)=xc;
     ycv(k)=yc; 
  else
     disp(sprintf('error openting file %s',fname))
  end

end

figure(3); clf;
subplot(2,1,1)
hold on;
title(pname)
subplot(2,1,1)
for i=1:8
  plot(timev,wlinev(i,:),[ccol(i)]);
end
ylabel('vorticity')
ax=axis;
axis([ax(1),ax(2),-1,ax(4)]);
hold off;

subplot(2,1,2)
hold on;
plot(timev,xcv-xcv(1),'b')
plot(timev,ycv-ycv(1),'r')
xlabel('time')
ylabel('x (blue), y (red)')
hold off;
orient tall
print -dpsc vxline.ps
print -djpeg vxline.jpg
return


