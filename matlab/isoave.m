
%name='/scratch1/taylorm/iso_256A0040.000'
%nx=256; delx_over_eta=1.73; epsilon=95.0;

%name='/scratch1/taylorm/iso12_250A0022.000'
%nx=250; delx_over_eta=.80; epsilon=3.9;

%name='../src/test0000.5000'
%nx=256; delx_over_eta=.7; epsilon=3.0;

name='../src/test0001.0000'
nx=256; delx_over_eta=2.00; epsilon=4.31;

%name='../src/test0001.5000'
%nx=256; delx_over_eta=3.46; epsilon=4.93;


fid=fopen([name,'.isostr'],'r','l');
eta = 1/(nx*delx_over_eta);
 
cdir=[ 'k','k','k' ];  % x,y,zz
cdir=[cdir, 'g','g','g','g','g','g'];  % face diagonals
cdir=[cdir, 'r','r','r','r'];      % body diagonals
cdir=[cdir, 'b','b','b','b','b','b','b','b','b','b','b','b'];      % 12 (1,2,0) directions
cdir=[cdir, 'y','y','y','y','y','y','y','y','y','y','y','y'];      % 12 (1,1,2) directions
cdir=[cdir, 'c','c','c','c','c','c','c','c','c','c','c','c'];      % 12 (1,2,2) directions



% get the weights:
w=textread('../src/voronoi/isoave.weights','%f');
% take every other weight
w=2*w(1:2:length(w));

%w=w./w/length(w);  % equally weighted


msize=4;   % marker size
xmax=1000;  % maximum x axis

ndelta=fread(fid,1,'float64');
ndir  =fread(fid,1,'float64');
nlon  =fread(fid,1,'float64');
ntran =fread(fid,1,'float64');
nnew1 =fread(fid,1,'float64');
nnew2 =fread(fid,1,'float64');

r_val=fread(fid,[ndelta,ndir],'float64');
r_val=r_val*delx_over_eta;            % convert to units of r/eta:
xx=(1:1.0:(nx./2.5))*delx_over_eta;   % units of r/eta
xx_box = xx/delx_over_eta/nx;         % in code units (box length)


if (nlon>=2) 
   D_ll=fread(fid,[ndelta,ndir],'float64');
   D_lll=fread(fid,[ndelta,ndir],'float64');
end
if (nlon>=4) 
   SP_lll=fread(fid,[ndelta,ndir],'float64');
   SN_lll=fread(fid,[ndelta,ndir],'float64');
end

if (ntran>=2) 
   D1_tt=fread(fid,[ndelta,ndir],'float64');
   D2_tt=fread(fid,[ndelta,ndir],'float64');
   D1_ltt=fread(fid,[ndelta,ndir],'float64');
   D2_ltt=fread(fid,[ndelta,ndir],'float64');
end
if (ntran>=4) 
   SP1_ltt=fread(fid,[ndelta,ndir],'float64');
   SP2_ltt=fread(fid,[ndelta,ndir],'float64');
   SN1_ltt=fread(fid,[ndelta,ndir],'float64');
   SN2_ltt=fread(fid,[ndelta,ndir],'float64');
end

% D_lll=SP_lll-SN_lll;
% but to make SP and SN have the same sign as D, multipl
% by -1 so:  D=SN-SP
SP_lll=-SP_lll;
SN_lll=-SN_lll;
SP1_ltt=-SP1_ltt;
SN1_ltt=-SN1_ltt;
SP2_ltt=-SP2_ltt;
SN2_ltt=-SN2_ltt;


disp('1 = Scaling laws for total structure function');
disp('2 = Scaling laws and also plot D+/-');
disp('4 = 2nd and 3rd order isotropy check');
disp('5 = 2nd and 3rd order isotropy check, x,y,z directions only');
in=input('Enter choice: ');

plot_posneg=0;
if (in==4)
   klaws=0;
   check_isotropy=1;
elseif (in==5)
   ndir=3;
   w=w/sum(w(1:ndir));
   klaws=0;
   check_isotropy=1;
else
   klaws=1;
   check_isotropy=0;
   if (in==2) 
     plot_posneg=1;
   end
end







if (klaws) 

%
%  the 4/5 law
%
figure(1)
yyave=0*xx;
yyave_sq=0*xx;

pave=yyave;
nave=yyave;

for i=1:ndir
  x=r_val(:,i);
  x_box = x/delx_over_eta/nx;  % in code units (box length)
  y=-D_lll(:,i)./(x_box*epsilon);

  semilogx(x,y,['.',cdir(i)],'MarkerSize',msize);   hold on;
  yy = spline(x,y,xx);
  %plot(xx,yy,[cdir(i)],'LineWidth',.2);
  
  yyave=yyave+w(i)*yy;
  yyave_sq=yyave_sq + w(i)*yy.^2;
  
  % positive and negative parts:
  y=-SP_lll(:,i)./(x_box*epsilon);
  pave=pave+w(i)*spline(x,y,xx);
  y=-SN_lll(:,i)./(x_box*epsilon);
  nave=nave+w(i)*spline(x,y,xx);
  
end
yyave_sq=sqrt(yyave_sq)/sqrt(ndir);

semilogx(xx,yyave,'k','LineWidth',1.0);
%errorbar(xx,yyave,yyave_sq);
title('D_{lll} / r\epsilon   (4/5 law)       blue: D-/D+');
xlabel('r/\eta');
x=1:xmax; semilogx(x,.8*x./x,'k');
if (plot_posneg)
  grid
  semilogx(xx,pave)
  semilogx(xx,nave)
end
ax=axis;  axis([1,xmax,ax(3),ax(4)]);
hold off;
print -dpsc 45.ps
print -djpeg 45.jpg





%
%  the 4/15 law
%
figure(2)
yyave1=0*xx;
yyave2=0*xx;
pave=yyave1;
nave=yyave1;
for i=1:ndir
  x=r_val(:,i);
  x_box = x/delx_over_eta/nx;  % in code units (box length)
  y1=-D1_ltt(:,i)./(x_box*epsilon);
  y2=-D2_ltt(:,i)./(x_box*epsilon);

  semilogx(x,y1,['.',cdir(i)],'MarkerSize',msize); hold on;
  semilogx(x,y2,['.',cdir(i)],'MarkerSize',msize);

  yyave1=yyave1+w(i)*spline(x,y1,xx);
  yyave2=yyave2+w(i)*spline(x,y2,xx);
  
  
  y1=-SP1_ltt(:,i)./(x_box*epsilon);
  pave=pave++w(i)*spline(x,y1,xx);
  y1=-SN1_ltt(:,i)./(x_box*epsilon);
  nave=nave++w(i)*spline(x,y1,xx);
  

end
semilogx(xx,yyave1,'r');
semilogx(xx,yyave2,'r');
semilogx(xx,.5*(yyave1+yyave2),'k','LineWidth',1.0);
title('D_{ltt} / r\epsilon  (4/15 law)');
xlabel('r/\eta');
if (plot_posneg)
  grid
  semilogx(xx,pave)
  semilogx(xx,nave)
end
x=1:xmax; semilogx(x,(4/15)*x./x,'k');
ax=axis;  axis([1,xmax,ax(3),ax(4)]);
hold off;
print -dpsc 415.ps
print -djpeg 415.jpg


%
%  the 4/3 law
%
yyave=0*xx;
figure(3)
yyave=0*xx;
pave=yyave;
nave=yyave;
for i=1:ndir
  x=r_val(:,i);
  x_box = x/delx_over_eta/nx;  % in code units (box length)
  y=-(D_lll(:,i) + D1_ltt(:,i) + D2_ltt(:,i))./(x_box*epsilon);

  semilogx(x,y,['.',cdir(i)],'MarkerSize',msize); hold on;
  yyave=yyave+w(i)*spline(x,y,xx);
  
  y=-(SP_lll(:,i) + SP1_ltt(:,i) + SP2_ltt(:,i))./(x_box*epsilon);
  pave=pave+w(i)*spline(x,y,xx);
  y=-(SN_lll(:,i) + SN1_ltt(:,i) + SN2_ltt(:,i))./(x_box*epsilon);
  nave=nave+w(i)*spline(x,y,xx);
  
end
if (plot_posneg)
  grid
  semilogx(xx,pave)
  semilogx(xx,nave)
end
semilogx(xx,yyave,'k','LineWidth',1.0);
title('4/3 law');
xlabel('r/\eta');
x=1:xmax; semilogx(x,(4/3)*x./x,'k');
xlabel('r/\eta');
ax=axis;  axis([1,xmax,ax(3),ax(4)]);
hold off;
print -dpsc 43.ps
print -djpeg 43.jpg

end





if (check_isotropy)

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
  
  semilogx(x,y1,['.',cdir(i)],'MarkerSize',msize); hold on
  semilogx(x,y2,['.',cdir(i)],'MarkerSize',msize);
  
  yyave1=yyave1 + w(i)*spline(x,y1,xx); 
  yyave2=yyave2 + w(i)*spline(x,y2,xx); 
  yyave=yyave + w(i)*spline(x,y,xx);  
  
end
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
title('D_{tt} (points)       angle ave(red)       D_{ll} + .5 r (D_{ll})'' (green)');
xlabel('r/\eta');
ylabel('scaled by r^{-2/3}')
xlabel('r/\eta');
ax=axis;  axis([1,xmax,ax(3),ax(4)]);
hold off;
print -dpsc isocheck.ps
print -djpeg isocheck.jpg



%
% Gotoh 3rd order style isotropy check
%
yyave=0*xx;
yyave1=0*xx;
yyave2=0*xx;
figure(2)
for i=1:ndir
  x=r_val(:,i);
  x_box = x/delx_over_eta/nx;  % in code units (box length)

  y  = -D_lll(:,i); 
  y1 = -D1_ltt(:,i)./x_box;
  y2 = -D2_ltt(:,i)./x_box;
  
  semilogx(x,y1,['.',cdir(i)],'MarkerSize',msize); hold on
  semilogx(x,y2,['.',cdir(i)],'MarkerSize',msize);
  
  yyave1=yyave1 + w(i)*spline(x,y1,xx); 
  yyave2=yyave2 + w(i)*spline(x,y2,xx); 
  yyave=yyave + w(i)*spline(x,y,xx);  
  
end
semilogx(xx,yyave1,'r');
semilogx(xx,yyave2,'r');


%
%  compute and plot: [ 1/6 d/dr r D_+lll ] /r
%
f = yyave.*xx_box/6;
l=length(f);
df = ( f(3:l)-f(1:l-2)) ./ (xx_box(3:l)-xx_box(1:l-2));
df = df./xx_box(2:l-1);
semilogx(xx(2:l-1),df,'g');
title('D_{ltt} (points)       angle ave(red)         (r D_{lll})''/6 (green)');
ylabel('scaled by 1/r')
xlabel('r/\eta');
ax=axis;  axis([1,xmax,ax(3),ax(4)]);
hold off;
print -dpsc isocheck3.ps
print -djpeg isocheck3.jpg
end
