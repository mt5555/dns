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



for i=1:ndir
  x = r_val(:,i);                   % units of box length
  x_plot=x*nx*delx_over_eta;  % units of r/eta

  y=-D_lll(:,i)./(x*epsilon);

  if (plot_points==1 & mod(i,3)==1) 
%     semilogx(x_plot,y,['.k-'],'MarkerSize',msize);   hold on;
  end     
  semilogx(x_plot,y,['k.'],'MarkerSize',msize);   hold on;
  yy = spline(x,y,xx);
  
  yyave=yyave+w(i)*yy;
  yyave_sq=yyave_sq + w(i)*yy.^2;
  
  y  = D_ll(:,i); 
  yyave1=yyave1 + w(i)*spline(x,y,xx)/(epsilon);
  
end
yyave_sq=sqrt(yyave_sq)/sqrt(ndir);
max(yyave)
plot(xx_plot,yyave,'k','LineWidth',2.0); hold on;
y45=yyave;

%title('D_{lll} / r\epsilon   (4/5 law) ');


yname=['   0';'0.25';'0.50';'0.75';' 1.0'];
ylabel('< (u(x+r)-u(x))^3 > / (\epsilon r)','FontSize',16);
xlabel('r/\eta','FontSize',16);

x=1:xmax; plot(x,(4/5)*x./x,'k');


axis([1,xmax,0,1])
hold off;
print('-deps',[bname,'_45.ps']);
print -djpeg 45.jpg



end