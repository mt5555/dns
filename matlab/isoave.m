mu=0;
ke=0;
nx=1;
delx_over_eta=1;
eta = 1/(nx*delx_over_eta);
ext='.isostr';


%name='/ccs/scratch/taylorm/decay/decay20480000.6034.new'
%nx=2048; delx_over_eta=-1; epsilon=-1;

%name='/home2/skurien/fractional_stats/sc1024A0001.4000.new'
%nx=1024; delx_over_eta=2.95; epsilon=3.57; teddy=1.05; % R_l=434

%name='/home/mataylo/codes/dns_data/sc1024A/sc1024A0002.0000.new'
%nx=1024; delx_over_eta=2.95; epsilon=3.57; teddy=1.05; % R_l=434

%name='/home2/skurien/helicity_data/helical_forced/hel256_hpi2/hel256_hpi2_0005.8000.new'
%nx=256; delx_over_eta=2.97; epsilon=2.72; teddy=1.24; % R_l=186

%name='/nh/nest/u/skurien/projects/helicity_data/helical_forced/hel512_hpi2/diag/skhel512a0009.0000.new'
%nx=512; delx_over_eta=2.5; epsilon=2.72; teddy=1.24; % R_l=250

%anisotropic structure functions have suffix .isostr4
%name='/auto/nest/u/skurien/dns/src/skhel512a0009.0000'
%nx=512; delx_over_eta=2.5; epsilon=2.72; teddy=1.24; % R_l=250
%ext='.isostr4';

%%anisotropic structure functions have suffix .isostr4
%name='/auto/nest/u/skurien/dns/src/skhel512a0006.0000'
%nx=512; delx_over_eta=2.5; epsilon=2.72; teddy=1.24; % R_l=250
%ext='.isostr4';

%anisotropic structure functions have suffix .isostr4
%name='/auto/nest/u/skurien/dns/src/skhel512a0005.0000'
%nx=512; delx_over_eta=2.5; epsilon=2.72; teddy=1.24; % R_l=250
%ext='.isostr4';

%potential vorticity correlation functions have suffix .bisostr
%name='/nh/u/skurien/projects/pv/data_analysis/lowforc/n210040.0000'
%nx=128; delx_over_eta=2.5; epsilon=2.72; teddy=1.24; % R_l=250
%Q_eps = 1;
ext='.bisostr';

name='/nh/u/skurien/projects/pv/data_analysis/lowforc/qg0170.0000'
nx=256; delx_over_eta=2.5; epsilon=2.72; teddy=1.24; % R_l=250
Q_eps = 1;
ext='.bisostr';

ndir_use=73;

%ndir_use=49;  disp('USING ONLY 49 DIRECTIONS')

disp('1 = Scaling laws for total structure function');
disp('2 = Scaling laws and also plot D+/-');
disp('4 = 2nd and 3rd order isotropy check');
disp('5 = 2nd and 3rd order isotropy check, x,y,z directions only');
disp('6 = Scaling laws for total structure function (for paper)');
disp('7 = plot 4th order functions');
disp('8 = plot helical structure functions');
disp('9 = Mixed (anisotropic) structure functions');
disp('10= Scaling of potential vorticity/velocity correlation');

in=input('Enter choice: ');


plot_posneg=0;
if (in==4)
  klaws=0;
  check_isotropy=1;
elseif (in==5)
  ndir_use=3;
  klaws=0;
  check_isotropy=1;
elseif (in==7)
   klaws=2;
   check_isotropy=0;
elseif (in==8)
   klaws=3;
   check_isotropy=0;
elseif (in==9)
   klaws=4;
   check_isotropy=0;
elseif (in==10)
   klaws=5
   check_isotropy=0;
else
   klaws=1;
   check_isotropy=0;
if (in==2) 
     plot_posneg=1;
     end
end

xx=(1:.5:(nx./2.5)) / nx;
%xx=(1:.5:100)/nx;

if (in==6) 
  compisoave_paper(name,ext,xx,ndir_use,klaws,plot_posneg,check_isotropy,1);
else
  [y45,y415,y43,eps,h_eps,y215]=compisoave(name,ext,xx,ndir_use,klaws,plot_posneg,check_isotropy,1);
end

return



  
xx_plot=xx*nx*delx_over_eta;
figure(10); hold off;
semilogx(xx_plot,y45,'k'); hold on;
axis([1,1000,0,1.5]);


yysum=0*y45;

for i=0:1
for j=0:1
for k=0:1
   ext = sprintf('.isostr%i%i%i',i,j,k);
if (klaws ==5)
   ext = sprintf('.bisostr');
end
  [xx,y45,y415,y43,y23,epsl]=compisoave(name,ext,ndir_use,klaws,plot_posneg,check_isotropy);
  yysum=yysum+y45*epsl/eps/8;
  figure(10);
  plot(xx,y45);
end
end
end
  figure(10);
  plot(xx,yysum,'r');  

hold off;







