
name='/scratch1/taylorm/iso_256A0040.000'
nx=256; delx_over_eta=1.73; epsilon=95.0;

name='/scratch1/taylorm/iso12_250A0022.000'
nx=250; delx_over_eta=.80; epsilon=3.9;

name='../src/temp0000.0000'
nx=64; delx_over_eta=1.0; epsilon=1.0;


fid=fopen([name,'.isostr'],'r','l');
eta = 1/(nx*delx_over_eta);
 
cdir=[ 'k','k','k' ];  % x,y,zz
cdir=[cdir, 'g','g','g','g','g','g'];  % face diagonals
cdir=[cdir, 'r','r','r','r'];      % body diagonals
cdir=[cdir, 'b','b','b','b','b','b','b','b','b','b','b','b'];      % 12 (1,2,0) directions
cdir=[cdir, 'y','y','y','y','y','y','y','y','y','y','y','y'];      % 12 (1,1,2) directions
cdir=[cdir, 'm','m','m','m','m','m','m','m','m','m','m','m'];      % 12 (1,1,2) directions



ndelta=fread(fid,1,'float64');
ndir  =fread(fid,1,'float64');
nlon  =fread(fid,1,'float64');
ntran =fread(fid,1,'float64');
nnew1 =fread(fid,1,'float64');
nnew2 =fread(fid,1,'float64');

r_val=fread(fid,[ndelta,ndir],'float64');
r_val=r_val*delx_over_eta;            % convert to units of r/eta:
xx=(1:.25:(nx./2.5))*delx_over_eta;   % units of r/eta
xx_box = xx/delx_over_eta/nx;         % in code units (box length)


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


%
% Gotoh style isotropy check
%
yyave=0*xx;
yyave1=0*xx;
yyave2=0*xx;
figure(1)
for i=1:ndir
  x=r_val(:,i);
  x_box = x/delx_over_eta/nx;  % in code units (box length)

  y  = D_ll(:,i); 
  y1 = D1_tt(:,i).*x_box.^(-2/3);
  y2 = D2_tt(:,i).*x_box.^(-2/3);
  
  semilogx(x,y1,['.',cdir(i)]); hold on
  semilogx(x,y2,['.',cdir(i)]);
  
  yyave1=yyave1 + spline(x,y1,xx); 
  yyave2=yyave2 + spline(x,y2,xx); 
  yyave=yyave + spline(x,y,xx);  
  
end
yyave=yyave/ndir; yyave1=yyave1/ndir; yyave2=yyave2/ndir;
semilogx(xx,yyave1,'r');
semilogx(xx,yyave2,'r');


%
%  compute and plot: (D_ll  + .5 r d/dr ( D_ll) )^(-2/3)
%
f = yyave;
l=length(f);
df = ( f(3:l)-f(1:l-2)) ./ (xx_box(3:l)-xx_box(1:l-2));
f2 = f(2:l-1) + .5*xx_box(2:l-1).*df;
f2 = f2 .* xx_box(2:l-1).^(-2/3);
semilogx(xx(2:l-1),f2,'g');
title('D_{tt} r^{-2/3} (blue)       angle ave(red)       F(D_{ll}) (green)');
hold off;
print -dpsc isocheck.ps



%
%  the 4/5 law
%
yyave=0*xx;
figure(2)
yyave=0*xx;
for i=1:ndir
  x=r_val(:,i);
  x_box = x/delx_over_eta/nx;  % in code units (box length)
  y=-D_lll(:,i)./(x_box*epsilon);

  semilogx(x,y,['.',cdir(i)]);   hold on;
  yyave=yyave+spline(x,y,xx);
end
yyave=yyave/ndir;
semilogx(xx,yyave,'r');
title('D_{lll} / r\epsilon   (4/5 law)');
semilogx(xx,.8*xx./xx,':');
hold off;



%
%  the 4/15 law
%
figure(3)
yyave1=0*xx;
yyave2=0*xx;
for i=1:ndir
  x=r_val(:,i);
  x_box = x/delx_over_eta/nx;  % in code units (box length)
  y1=-D1_ltt(:,i)./(x_box*epsilon);
  y2=-D2_ltt(:,i)./(x_box*epsilon);

  semilogx(x,y1,['.',cdir(i)]); hold on;
  semilogx(x,y2,['.',cdir(i)]);

  yyave1=yyave1+spline(x,y1,xx);
  yyave2=yyave2+spline(x,y2,xx);
end
yyave1=yyave1/ndir;
yyave2=yyave2/ndir;
semilogx(xx,yyave1,'r');
semilogx(xx,yyave2,'r');
title('D_{ltt} / r\epsilon  (4/15 law)');
semilogx(xx,(4/15)*xx./xx,':');
hold off;



%
%  the 4/3 law
%
yyave=0*xx;
figure(4)
yyave=0*xx;
for i=1:ndir
  x=r_val(:,i);
  x_box = x/delx_over_eta/nx;  % in code units (box length)
  y=-(D_lll(:,i) + D1_ltt(:,i) + D2_ltt(:,i))./(x_box*epsilon);

  semilogx(x,y,['.',cdir(i)]); hold on;
  yyave=yyave+spline(x,y,xx);
end
yyave=yyave/ndir;
semilogx(xx,yyave,'r');
title('4/3 law');
semilogx(xx,(4/3)*xx./xx,':');
hold off;



