function [y45,y415,y43,epsilon]=compisoave(name,ext,xx,ndir_use,klaws,plot_posneg,check_isotropy,plot_points)


l=findstr('/',name);
l=l(length(l));
bname=name(l+1:length(name));
if (plot_posneg) bname=[bname,'_s']; end;
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

[nx,ndelta,ndir,r_val,ke,epsilon,mu,...
    D_ll,D_lll,D1_tt,D2_tt,D1_ltt,D2_ltt,...
    SP_lll,SN_lll,SP1_ltt,SP2_ltt,SN1_ltt,SN2_ltt,H_ltt,H_tt] ...
     = readisostr( [name,ext] );


eta = (mu^3 / epsilon)^.25;
delx_over_eta=(1/nx)/eta;


%
% use only 49 directions:
if (ndir_use>0) ndir=ndir_use; end;

if (ndir==3)
  w=ones([1,3])/3;
else
  equalw=1;
  if (equalw) 
    % put this in to use equally weighted:
    w=ones([1,ndir])/ndir;
  else
    % get the weights:
    wname=sprintf('../src/voronoi/isoave.weights %i',ndir);
    w=textread(wname,'                  %f');
    % take every other weight
    w=2*w(1:2:length(w));
  end
end
if (abs(1-sum(w))>1e-7) 
  disp('error: weights do not sum to 1')
  return;
end




%
% xx is given in units of the box length
% but r_val is given in units of delx.  convert to box length units
%

r_val=r_val/nx;
xx_plot = xx*nx*delx_over_eta;        % units of r/eta

lambda=sqrt(10*ke*mu/epsilon);       % single direction lambda
R_lambda = lambda*sqrt(2*ke/3)/mu;


disp(sprintf('KE:      %f  2pi units: %f',ke,ke*4*pi*pi));
disp(sprintf('epsilon: %f  2pi units: %f',epsilon,epsilon*4*pi*pi));
disp(sprintf('mu:      %f  2pi units: %f',mu,mu*4*pi*pi));
disp(sprintf('eta:     %f  2pi units: %f',eta,eta*2*pi));
disp(sprintf('lambda:  %f  2pi units: %f',lambda,lambda*2*pi));
disp(sprintf('delx/eta %f',delx_over_eta));
disp(sprintf('R_l:     %f',R_lambda));
disp(sprintf('ndir:     %f',ndir));
disp(' ')






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
  x = r_val(:,i);                   % units of box length
  x_plot=x*nx*delx_over_eta;  % units of r/eta

  y=-D_lll(:,i)./(x*epsilon);

  if (plot_points==1) 
     semilogx(x_plot,y,['.',cdir(i)],'MarkerSize',msize);   hold on;
  end     
  yy = spline(x,y,xx);
  
  yyave=yyave+w(i)*yy;
  yyave_sq=yyave_sq + w(i)*yy.^2;
  
  % positive and negative parts:
  y=-SP_lll(:,i)./(x*epsilon);
  pave=pave+w(i)*spline(x,y,xx);
  y=-SN_lll(:,i)./(x*epsilon);
  nave=nave+w(i)*spline(x,y,xx);

  y  = D_ll(:,i); 
  yyave1=yyave1 + w(i)*spline(x,y,xx)/(epsilon);
  
end
yyave_sq=sqrt(yyave_sq)/sqrt(ndir);
max(yyave)
plot(xx_plot,yyave,'k','LineWidth',1.0); hold on;
y45=yyave;

title('D_{lll} / r\epsilon   (4/5 law) ');
ylabel(pname);
xlabel('r/\eta');
x=1:xmax; plot(x,(4/5)*x./x,'k');
if (plot_posneg)
  grid
  plot(xx_plot,pave)
  plot(xx_plot,nave)
end


ax=axis;  axis([1,xmax,ax(3),ax(4)]);
hold off;
if (plot_points==1) 
print('-dpsc',[bname,'_45.ps']);
print -djpeg 45.jpg
end

if (0) 
figure(4)
hold off; clf;
%curve fit (xx,pave)  and (xx,nave)
%  
%  pave, nave have been divied by x
%  compute: 6 mu (1/r) d/dr yyave1
l=length(yyave1);
df = ( yyave1(3:l)-yyave1(1:l-2)) ./ (xx(3:l)-xx(1:l-2));
df = 6*mu*df./xx(3:l);

x=xx_plot(3:l);
pave=pave(3:l);
nave=nave(3:l);
semilogx(x,pave+.5*df,'k',x,nave+.5*df,'k',x,pave,x,nave);
grid;

%Bm = lsqcurvefit('fun3',[0],xx,nave)
%Bp = lsqcurvefit('fun3',[0],xx,pave)
%bm = fun3([Bm(1)],xx);
%bp = fun3([Bp(1)],xx);
%plot(xx,bm,xx,bp)
end



%
%  the 4/15 law
%
figure(2)
yyave1=0*xx;
yyave2=0*xx;
pave=yyave1;
nave=yyave1;
for i=1:ndir
  x = r_val(:,i);                       % units of box length
  x_plot=x*nx*delx_over_eta;  % units of r/eta

  y1=-D1_ltt(:,i)./(x*epsilon);
  y2=-D2_ltt(:,i)./(x*epsilon);

  if (plot_points==1) 
  semilogx(x_plot,y1,['.',cdir(i)],'MarkerSize',msize); hold on;
  semilogx(x_plot,y2,['.',cdir(i)],'MarkerSize',msize);
  end
  
  yyave1=yyave1+w(i)*spline(x,y1,xx);
  yyave2=yyave2+w(i)*spline(x,y2,xx);
  
  y1=-.5*(SP2_ltt(:,i)+SP1_ltt(:,i))./(x*epsilon);
  pave=pave++w(i)*spline(x,y1,xx);
  y1=-.5*(SN2_ltt(:,i)+SN1_ltt(:,i))./(x*epsilon);
  nave=nave++w(i)*spline(x,y1,xx);
  

end
%semilogx(xx_plot,yyave1,'r'); hold on;
%semilogx(xx_plot,yyave2,'r');
semilogx(xx_plot,.5*(yyave1+yyave2),'k','LineWidth',1.0);
title('D_{ltt} / r\epsilon  (4/15 law)');
y415=.5*(yyave1+yyave2);
ylabel(pname);
xlabel('r/\eta');
if (plot_posneg)
  grid
  semilogx(xx_plot,pave)
  semilogx(xx_plot,nave)
end
x=1:xmax; semilogx(x,(4/15)*x./x,'k');
ax=axis;  axis([1,xmax,ax(3),ax(4)]);
hold off;
if (plot_points==1) 
print('-dpsc',[bname,'_415.ps']);
print -djpeg 415.jpg
end

%
%  the 4/3 law
%
yyave=0*xx;
figure(3)
yyave=0*xx;
pave=yyave;
nave=yyave;
for i=1:ndir
  x = r_val(:,i);                       % units of box length
  x_plot=x*nx*delx_over_eta;  % units of r/eta

  y=-(D_lll(:,i) + D1_ltt(:,i) + D2_ltt(:,i))./(x*epsilon);

  if (plot_points==1) 
    semilogx(x_plot,y,['.',cdir(i)],'MarkerSize',msize); hold on;
  end      
  yyave=yyave+w(i)*spline(x,y,xx);
  
  y=-(SP_lll(:,i) + SP1_ltt(:,i) + SP2_ltt(:,i))./(x*epsilon);
  pave=pave+w(i)*spline(x,y,xx);
  y=-(SN_lll(:,i) + SN1_ltt(:,i) + SN2_ltt(:,i))./(x*epsilon);
  nave=nave+w(i)*spline(x,y,xx);
  
end
semilogx(xx_plot,yyave,'k','LineWidth',1.0); hold on;
if (plot_posneg)
  grid
  semilogx(xx_plot,pave)
  semilogx(xx_plot,nave)
end

y43=yyave1;
title('4/3 law');
ylabel(pname);
xlabel('r/\eta');
x=1:xmax; semilogx(x,(4/3)*x./x,'k');
xlabel('r/\eta');
ax=axis;  axis([1,xmax,ax(3),ax(4)]);
hold off;
if (plot_points==1) 
print('-dpsc',[bname,'_43.ps']);
print -djpeg 43.jpg
end


%
%  the 2/15 law
%
if (0)
yyave=0*xx;
figure(4)
yyave=0*xx;
pave=yyave;
nave=yyave;
for i=1:ndir
  x = r_val(:,i);                       % units of box length
  x_plot=x*nx*delx_over_eta;            % units of r/eta

  y=H_ltt(:,i)./(epsilon*(x.^2));

  semilogx(x_plot,y,['o-',cdir(i)],'MarkerSize',msize); hold on;
  yyave=yyave+w(i)*spline(x,y,xx);

end

semilogx(xx_plot,yyave,'k','LineWidth',2.5);
y215=yyave;
title('2/15 law');
ylabel(pname);
xlabel('r/\eta');
x=1:xmax; semilogx(x,(2/15)*x./x,'k');
xlabel('r/\eta');
ax=axis;  axis([1,xmax,ax(3),ax(4)]);
hold off;
if (plot_points==1) 
print('-dpsc',[bname,'_215.ps']);
print -djpeg 215.jpg
end
end




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
  x = r_val(:,i);                       % units of box length
  x_plot=x*nx*delx_over_eta;  % units of r/eta

  y  = D_ll(:,i); 
  y1 = D1_tt(:,i).*x.^(-2/3);
  y2 = D2_tt(:,i).*x.^(-2/3);
  
  semilogx(x_plot,y1,['.',cdir(i)],'MarkerSize',msize); hold on
  semilogx(x_plot,y2,['.',cdir(i)],'MarkerSize',msize);
  
  yyave1=yyave1 + w(i)*spline(x,y1,xx); 
  yyave2=yyave2 + w(i)*spline(x,y2,xx); 
  yyave=yyave + w(i)*spline(x,y,xx);  
  
end
semilogx(xx_plot,.5*(yyave1+yyave2),'r');


%
%  compute and plot: (D_ll  + .5 r d/dr ( D_ll) )^(-2/3)
%
f = yyave;
l=length(f);
df = ( f(3:l)-f(1:l-2)) ./ (xx(3:l)-xx(1:l-2));
f2 = f(2:l-1) + .5*xx(2:l-1).*df;
f2 = f2 .* xx(2:l-1).^(-2/3);
semilogx(xx_plot(2:l-1),f2,'g');
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
  x = r_val(:,i);                       % units of box length
  x_plot=x*nx*delx_over_eta;  % units of r/eta

  y  = -D_lll(:,i); 
  y1 = -D1_ltt(:,i)./x;
  y2 = -D2_ltt(:,i)./x;
  
  semilogx(x_plot,y1,['.',cdir(i)],'MarkerSize',msize); hold on
  semilogx(x_plot,y2,['.',cdir(i)],'MarkerSize',msize);
  
  yyave1=yyave1 + w(i)*spline(x,y1,xx); 
  yyave2=yyave2 + w(i)*spline(x,y2,xx); 
  yyave=yyave + w(i)*spline(x,y,xx);  
  
end
semilogx(xx_plot,.5*(yyave1+yyave2),'r');


%
%  compute and plot: [ 1/6 d/dr r D_+lll ] /r
%
f = yyave.*xx/6;
l=length(f);
df = ( f(3:l)-f(1:l-2)) ./ (xx(3:l)-xx(1:l-2));
df = df./xx(2:l-1);
semilogx(xx_plot(2:l-1),df,'g');
title('D_{ltt} (points)       angle ave(red)         (r D_{lll})''/6 (green)');
ylabel(pname);
xlabel('r/\eta');
ax=axis;  axis([1,xmax,ax(3),ax(4)]);
hold off;
print('-dpsc',[bname,'_isocheck3.ps']);
print -djpeg isocheck3.jpg
end
