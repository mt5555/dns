function [xx_out,y45,y415,y43,epsilon]=compisoave(name,ext,ndir_use,klaws,plot_posneg,check_isotropy)

fid=fopen([name,ext],'r','l');

l=findstr('/',name);
l=l(length(l));
bname=name(l+1:length(name));
l=findstr('_',bname);
pname=[bname,ext];
pname(l)='-';

cdir=[ 'k','k','k' ];  % x,y,zz
cdir=[cdir, 'g','g','g','g','g','g'];  % face diagonals
cdir=[cdir, 'r','r','r','r'];      % body diagonals
cdir=[cdir, 'b','b','b','b','b','b','b','b','b','b','b','b'];      % 12 (1,2,0) directions
cdir=[cdir, 'b','b','b','b','b','b','b','b','b','b','b','b'];      % 12 (1,1,2) directions
cdir=[cdir, 'b','b','b','b','b','b','b','b','b','b','b','b'];      % 12 (1,2,2) directions
cdir=[cdir, 'y','y','y','y','y','y','y','y','y','y','y','y'];      % 12 (1,3,0) directions
cdir=[cdir, 'y','y','y','y','y','y','y','y','y','y','y','y'];      % 12 (1,1,3) directions





msize=4;   % marker size
xmax=1000;  % maximum x axis

ndelta=fread(fid,1,'float64');
ndir  =fread(fid,1,'float64');
nlon  =fread(fid,1,'float64');
ntran =fread(fid,1,'float64');
nscalars =fread(fid,1,'float64');
nnew2 =fread(fid,1,'float64');


r_val=fread(fid,[ndelta,ndir],'float64');
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


if (nscalars==7) 
  time=fread(fid,1,'float64');
  nx=fread(fid,1,'float64');
  ny=fread(fid,1,'float64');
  nz=fread(fid,1,'float64');
  mu=fread(fid,1,'float64');
  ke=fread(fid,1,'float64');
  epsilon=fread(fid,1,'float64');
  eta = (mu^3 / epsilon)^.25;
  delx_over_eta=(1/nx)/eta;
end


%
% use only 49 directions:
if (ndir_use>0) ndir=ndir_use; end;

if (ndir==3)
  w=ones([1,3])/3;
else
  % get the weights:
  wname=sprintf('../src/voronoi/isoave.weights%i',ndir);
  w=textread(wname,'%f');
  % take every other weight
  w=2*w(1:2:length(w));

  % put this in to use equally weighted:
  %w=w./w/length(w);  % equally weighted
end
if (abs(1-sum(w))>1e-7) 
  disp('error: weights do not sum to 1')
  return;
end





r_val=r_val*delx_over_eta;            % convert to units of r/eta:
xx=(1:.5:(nx./2.5))*delx_over_eta;   % units of r/eta
xx_box = xx/delx_over_eta/nx;         % in code units (box length)

lambda=sqrt(10*ke*mu/epsilon);       % single direction lambda
R_lambda = lambda*sqrt(2*ke/3)/mu;


disp(sprintf('KE:      %f',ke));
disp(sprintf('epsilon: %f',epsilon));
disp(sprintf('mu:      %f',mu));
disp(sprintf('eta:     %f',eta));
disp(sprintf('delx/eta %f',delx_over_eta));
disp(sprintf('lambda:  %f',lambda));
disp(sprintf('R_l:     %f',R_lambda));
disp(' ')



% D_lll=SP_lll-SN_lll;
% but to make SP and SN have the same sign as D, multipl
% by -1 so:  D=SN-SP
SP_lll=-SP_lll;
SN_lll=-SN_lll;
SP1_ltt=-SP1_ltt;
SN1_ltt=-SN1_ltt;
SP2_ltt=-SP2_ltt;
SN2_ltt=-SN2_ltt;








if (klaws) 

%
%  the 4/5 law
%
figure(1)
yyave=0*xx;
yyave_sq=0*xx;
yyave1=yyave;

y45=yyave;
y415=yyave;
y43=yyave;

pave=yyave;
nave=yyave;

for i=1:ndir
  x=r_val(:,i);
  x_box = x/delx_over_eta/nx;  % in code units (box length)
  y=-D_lll(:,i)./(x_box*epsilon*4/5);

  semilogx(x,y,['.',cdir(i)],'MarkerSize',msize);   hold on;
  %loglog(x,y,['.',cdir(i)],'MarkerSize',msize);   hold on;
  yy = spline(x,y,xx);
  %plot(xx,yy,[cdir(i)],'LineWidth',.2);
  
  yyave=yyave+w(i)*yy;
  yyave_sq=yyave_sq + w(i)*yy.^2;
  
  % positive and negative parts:
  y=-SP_lll(:,i)./(x_box*epsilon*4/5);
  pave=pave+w(i)*spline(x,y,xx);
  y=-SN_lll(:,i)./(x_box*epsilon*4/5);
  nave=nave+w(i)*spline(x,y,xx);

  y  = D_ll(:,i); 
  yyave1=yyave1 + w(i)*spline(x,y,xx)/(epsilon*4/5);
  
end
yyave_sq=sqrt(yyave_sq)/sqrt(ndir);
max(yyave*4/5)
plot(xx,yyave,'k','LineWidth',1.0);
y45=yyave;
xx_out=xx;
%errorbar(xx,yyave,yyave_sq);
title('D_{lll} / 0.8 r\epsilon   (4/5 law) ');
ylabel(pname);
xlabel('r/\eta');
x=1:xmax; plot(x,1.0*x./x,'k');
if (plot_posneg)
  grid
  plot(xx,pave)
  plot(xx,nave)
end


ax=axis;  axis([1,xmax,-.5,2.0]);
hold off;
print('-dpsc',[bname,'_45.ps']);
print -djpeg 45.jpg

figure(4)
%curve fit (xx,pave)  and (xx,nave)
%  
%  pave, nave have been divied by x
%  compute: 6 mu (1/r) d/dr yyave1
l=length(yyave1);
df = ( yyave1(3:l)-yyave1(1:l-2)) ./ (xx_box(3:l)-xx_box(1:l-2));
df = 6*mu*df./xx_box(3:l);

xx=xx(3:l);
pave=pave(3:l);
nave=nave(3:l);
semilogx(xx,pave+.5*df,'k',xx,nave+.5*df,'k',xx,pave,xx,nave);
grid;

%Bm = lsqcurvefit('fun3',[0],xx,nave)
%Bp = lsqcurvefit('fun3',[0],xx,pave)
%bm = fun3([Bm(1)],xx);
%bp = fun3([Bp(1)],xx);
%plot(xx,bm,xx,bp)




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
  y1=-D1_ltt(:,i)./(x_box*epsilon*4/15);
  y2=-D2_ltt(:,i)./(x_box*epsilon*4/15);

  semilogx(x,y1,['.',cdir(i)],'MarkerSize',msize); hold on;
  semilogx(x,y2,['.',cdir(i)],'MarkerSize',msize);

  yyave1=yyave1+w(i)*spline(x,y1,xx);
  yyave2=yyave2+w(i)*spline(x,y2,xx);
  
  
  y1=-SP1_ltt(:,i)./(x_box*epsilon*4/15);
  pave=pave++w(i)*spline(x,y1,xx);
  y1=-SN1_ltt(:,i)./(x_box*epsilon*4/15);
  nave=nave++w(i)*spline(x,y1,xx);
  

end
semilogx(xx,yyave1,'r');
semilogx(xx,yyave2,'r');
semilogx(xx,.5*(yyave1+yyave2),'k','LineWidth',1.0);
title('D_{ltt} / 4/15 r\epsilon  (4/15 law)');
y415=.5*(yyave1+yyave2);
ylabel(pname);
xlabel('r/\eta');
if (plot_posneg)
  grid
  semilogx(xx,pave)
  semilogx(xx,nave)
end
x=1:xmax; semilogx(x,x./x,'k');
ax=axis;  axis([1,xmax,-.5,2.0]);
hold off;
print('-dpsc',[bname,'_415.ps']);
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
  y=-(D_lll(:,i) + D1_ltt(:,i) + D2_ltt(:,i))./(x_box*epsilon*4/3);

  semilogx(x,y,['.',cdir(i)],'MarkerSize',msize); hold on;
  yyave=yyave+w(i)*spline(x,y,xx);
  
  y=-(SP_lll(:,i) + SP1_ltt(:,i) + SP2_ltt(:,i))./(x_box*epsilon*4/3);
  pave=pave+w(i)*spline(x,y,xx);
  y=-(SN_lll(:,i) + SN1_ltt(:,i) + SN2_ltt(:,i))./(x_box*epsilon*4/3);
  nave=nave+w(i)*spline(x,y,xx);
  
end
if (plot_posneg)
  grid
  semilogx(xx,pave)
  semilogx(xx,nave)
end
semilogx(xx,yyave,'k','LineWidth',1.0);
y43=yyave1;
title('4/3 law');
ylabel(pname);
xlabel('r/\eta');
x=1:xmax; semilogx(x,x./x,'k');
xlabel('r/\eta');
ax=axis;  axis([1,xmax,-.5,2.0]);
hold off;
print('-dpsc',[bname,'_43.ps']);
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
ylabel(pname);
xlabel('r/\eta');
ax=axis;  axis([1,xmax,ax(3),ax(4)]);
hold off;
print('-dpsc',[bname,'_isocheck2.ps']);
print -djpeg isocheck2.jpg



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
ylabel(pname);
xlabel('r/\eta');
ax=axis;  axis([1,xmax,ax(3),ax(4)]);
hold off;
print('-dpsc',[bname,'_isocheck3.ps']);
print -djpeg isocheck3.jpg
end
