%
%########################################################################
%#  plot passive scalar
%########################################################################
%


name = '/scratch2/taylorm/tmix256C/tmix256C'
time=1.0984;
times=sprintf('%.5f',time+10000);
times=times(2:length(times)-1);

npassive=10;
for np=1:npassive
  ext=sprintf('%3i',np+100);
  ext=ext(2:length(ext));
  ext=['.s',ext];
  
  s=findstr(name,'/');
  s=s(length(s));
  shortname=name(s+1:length(name));

   fname=[name,times,ext]
   [x,y,z,s,time]=getfield(fname);
   smax=max(max(s));
   [mx,slice1]=max(smax);
   smax=min(min(s));
   [mn,slice2]=min(smax);


   subplot(npassive/2,2,np)
   splot=squeeze(s(:,:,slice1));
   pcolor(x,y,splot')
   
   stitle=sprintf('%s    time=%.2f  max=%f',shortname,time,mx)
   if (np==1) title(stitle); end;
   axis equal
   axis([0,max(x),0,max(y)]);
   shading interp
end
   orient tall
   print('-djpeg','-r200',['p',times,'.jpg']); 

