% fracisoave  (name,ext,xx,ndir_use)
%
% compute angle average and plot a single snapshot
%
%


name='/scratch2/taylorm/sc1024A/sc1024A0001.9000'
ext='.isostrf';
nx=1024; delx_over_eta=2.95; epsilon=3.57; teddy=1.05; % R_l=434



%  index=1 .. 9   corresponds to u_l^(.1*p)
p=5; 
frac_power = .1*p;
plot_points=1;
ndir_use=0;

xx=(1:.5:(nx./2.5)) / nx;





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
xmax=3000;  % maximum x axis
iso_check_dir=2;  % direction to use for single direction iso_check


[nx,ndelta,ndir,r_val,ke,epsilon,mu,...
    D_ll,D_lll,D1_tt,D2_tt,D1_ltt,D2_ltt,...
    SP_lll,SN_lll,SP1_ltt,SP2_ltt,SN1_ltt,SN2_ltt,H_ltt,H_tt,D_lltt,Dl,Dt,...
    h_epsilon] ...
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
disp(sprintf('h_epsilon: %f  2pi units: %f',h_epsilon,h_epsilon*2*pi));
disp(sprintf('mu:      %f  2pi units: %f',mu,mu*4*pi*pi));
disp(sprintf('eta:     %f  2pi units: %f',eta,eta*2*pi));
disp(sprintf('lambda:  %f  2pi units: %f',lambda,lambda*2*pi));
disp(sprintf('delx/eta %f',delx_over_eta));
disp(sprintf('R_l:     %f',R_lambda));
disp(sprintf('ndir:     %f',ndir));
disp(' ')



figure(1); subplot(1,1,1);
yyave_t=0*xx;
yyave_l=0*xx;


ndir_vec=1:ndir;



figure(1); clf; subplot(1,1,1)
for i=ndir_vec
  x = r_val(:,i);                   % units of box length
  x_plot=x*nx*delx_over_eta;  % units of r/eta
  
  yl=Dl(:,i,p);
  
  if (plot_points==1) 
    semilogx(x_plot,yl,['.',cdir(i)],'MarkerSize',msize);   hold on; 
  end     
  yyl = spline(x,yl,xx);
  
  yyave_l=yyave_l+w(i)*yyl;
  
end


title(sprintf('D_l^{%.1f}',frac_power));
ylabel(pname);
xlabel('r/\eta');
plot(xx_plot,yyave_l,'k','LineWidth',1.0); hold on;
ax=axis;  axis([1,xmax,ax(3),ax(4)]);

hold off;




figure(2); clf; subplot(1,1,1)
for i=ndir_vec
  x = r_val(:,i);                   % units of box length
  x_plot=x*nx*delx_over_eta;  % units of r/eta

  yt=Dt(:,i,p);

  if (plot_points==1) 
     semilogx(x_plot,yt,['.',cdir(i)],'MarkerSize',msize);   hold on;
  end     
  yyt = spline(x,yt,xx);
  
  yyave_t=yyave_t+w(i)*yyt;
  
end


title(sprintf('D_t^{%.1f}',frac_power));
ylabel(pname);
xlabel('r/\eta');
plot(xx_plot,yyave_t,'k','LineWidth',1.0); hold on;
ax=axis;  axis([1,xmax,ax(3),ax(4)]);

hold off;




