name='/ccs/scratch/taylorm/dns/iso12_512'
name='/home/mt/iso12_512'
nx=512; delx_over_eta=2.74; epsilon=3.89;  teddy=1.0;

ext='.isostr';
ext2='.isow2s2';

times=[7.0:.1:7.0];
xx=(1:.5:(nx./2.5)) / nx;
xx_plot=(1:.5:(nx./2.5)) *delx_over_eta;   % units of r/eta

t=times(1);
tstr=sprintf('%10.4f',t+10000);
fname=[name,tstr(2:10)];
disp([fname,ext]);


[nx,ndelta,ndir,r_val,ke,eps_l,mu,...
    D_ll,D_lll,D1_tt,D2_tt,D1_ltt,D2_ltt,...
    SP_lll,SN_lll,SP1_ltt,SP2_ltt,SN1_ltt,SN2_ltt,H_ltt,H_tt,D_lltt,Dl] ...
     = readisostr([fname,ext]);
[ndelta_2,ndir_2,r_val_2,w2w2,s2s2,w2s2] = readisow2s2([fname,ext2]);
    
    
eta_l = (mu^3 / eps_l)^.25;
delx_over_eta_l=(1/nx)/eta_l;
    
% convert to units of box length:
r_val=r_val/nx;
r_val2=r_val2/nx;

msize=4;   % marker size
xmax=1000;  % maximum x axis



%
%  angle average
%
figure(2)
w2w2ave=0*xx;
s2s2ave=0*xx;
w2s2ave=0*xx;
D_llttave=0*xx;
D_llllave=0*xx;

for i=1:ndir
  x = r_val(:,i);                 % units of box length
  x_plot=x*nx*delx_over_eta_l;    % units of r/eta

  %Dl=Dl(ndelta,ndir,p-1);  for <del_u**p>  p=2..10
  p=4;
  y1=Dl(:,i,p-1);
  semilogx(x_plot,y1,['.',cdir(i)],'MarkerSize',msize); hold on;
  D_llllave=D_llllave+w(i)*spline(x,y1,xx);
  
  y1=D_lltt(:,i);
  semilogx(x_plot,y1,['.',cdir(i)],'MarkerSize',msize); hold on;
  D_llttave=D_llttave+w(i)*spline(x,y1,xx);


  % this data also includes r=0
  x = r_val2(:,i);                       % units of box length
  x_plot2=x*nx*delx_over_eta_l;          % units of r/eta

  y1=w2w2(:,i);
  semilogx(x_plot,y1,['.',cdir(i)],'MarkerSize',msize); hold on;
  w2w2ave=w2w2ave+w(i)*spline(x,y1,xx);

  y1=s2s2(:,i);
  semilogx(x_plot,y1,['.',cdir(i)],'MarkerSize',msize); hold on;
  s2s2ave=s2s2ave+w(i)*spline(x,y1,xx);
  
  y1=w2s2(:,i);
  semilogx(x_plot,y1,['.',cdir(i)],'MarkerSize',msize); hold on;
  w2s2ave=w2s2ave+w(i)*spline(x,y1,xx);
end

semilogx(xx_plot,w2w2ave,'k','LineWidth',1.0);

title(' ');
ylabel(' ');
xlabel('r/\eta');
%x=1:xmax; semilogx(x,(4/15)*x./x,'k');
ax=axis;  axis([1,xmax,ax(3),ax(4)]);
hold off;

%print('-dpsc',[name,ext,'_w2s2.ps']);
%print -djpeg w2s2.jpg


