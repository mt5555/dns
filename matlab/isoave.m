
name='../src/temp0000.0000';
fid=fopen([name,'.isostr'],'r','l');


ndelta=fread(fid,1,'float64');
ndir  =fread(fid,1,'float64');
nlon  =fread(fid,1,'float64');
ntran =fread(fid,1,'float64');
nnew1 =fread(fid,1,'float64');
nnew2 =fread(fid,1,'float64');

r_val=fread(fid,[ndelta,ndir],'float64');
if (nlon==2) 
   D_ll=fread(fid,[ndelta,ndir],'float64');
   D_lll=fread(fid,[ndelta,ndir],'float64');
end
if (ntran==2) 
   D1_tt=fread(fid,[ndelta,ndir],'float64');
   D2_tt=fread(fid,[ndelta,ndir],'float64');
   D1_ltt=fread(fid,[ndelta,ndir],'float64');
   D2_ltt=fread(fid,[ndelta,ndir],'float64');
end

xx=1:.25:32;
yyave=0*xx;

figure(1)
for i=1:ndir
  x=r_val(:,i); 
  y=D_ll(:,i);
  semilogx( x,y,'.' )
  hold on;
  yy=spline(x,y,xx);
  %semilogx( xx,yy);
  yyave=yyave+yy;
end
yyave=yyave/ndir;
semilogx( xx,yyave,'r');
%axis([1 64 0 4]);
hold off;


figure(2)
yyave=0*xx;
for i=1:ndir
  x=r_val(:,i); 
  y=D_lll(:,i)./(x.^3);
  semilogx( x,y,'.' )
  hold on;
  yy=spline(x,y,xx);
  %semilogx( xx,yy);
  yyave=yyave+yy;
end
yyave=yyave/ndir;
semilogx( xx,yyave,'r');
%axis([1 64 -.01 .01]);
hold off;



