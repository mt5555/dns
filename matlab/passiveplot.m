%
%########################################################################
%#  plot passive scalar
%########################################################################
%


name = '/scratch2/taylorm/tmix256C/tmix256C'
ext='.s01';
times=[1.0509];

s=findstr(name,'/');
s=s(length(s));
shortname=name(s+1:length(name));



   ts = sprintf('%9.5f',10000+times);
   ts=ts(2:10);
   fname=[name,ts,ext]
   [x,y,z,s,time]=getfield(fname);
   smax=max(max(s.^2));
   [m1,slice1]=max(smax);
   [m1,slice2]=min(smax);


   subplot(2,1,1)
   splot=squeeze(s(:,:,slice1));
   pcolor(x,y,splot')

   stitle=sprintf('%s    time=%.2f  max=%f',shortname,time)
   title(stitle);
   axis equal
   axis([0,max(x),0,max(y)]);
   shading interp

   subplot(2,1,2)
   splot=squeeze(s(:,:,slice1));
   pcolor(x,y,splot')
   axis equal
   axis([0,max(x),0,max(y)]);
   shading interp
