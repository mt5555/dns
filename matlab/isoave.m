mu=0;
ke=0;
nx=1;
delx_over_eta=1;
eta = 1/(nx*delx_over_eta);
ext='.isostr';

%name='/scratch1/taylorm/iso12w512A0001.3847'
%nx=512; delx_over_eta=5.8615; epsilon=.2849;

%name='/scratch1/taylorm/iso12_500A0001.7723'
%nx=500; delx_over_eta=2.740; epsilon=3.5208;

%name='/scratch1/taylorm/iso12_250A0022.000'
%nx=250; delx_over_eta=.80; epsilon=3.9;


name='/ccs/scratch/taylorm/dns/iso12/iso12_5120002.7000'
nx=512; delx_over_eta=2.75; epsilon=3.95;  %R_l=249
%ext='.isostr001';

%name='/ccs/scratch/taylorm/sk/check256_0000.8000'
%nx=256; delx_over_eta=1.0; epsilon=1.0;

%name='/ccs/scratch/taylorm/decay/decay2048-1024.0000.6491'
%nx=1024; delx_over_eta=2.73*2; epsilon=.04;

%name='/ccs/scratch/taylorm/decay/decay20480000.6034.new'
%nx=2048; delx_over_eta=-1; epsilon=-1;

%name='/ccs/scratch/taylorm/dns/sc1024A/sc1024A0002.0000.new'
%nx=1024; delx_over_eta=2.95; epsilon=3.57; teddy=1.05; % R_l=434


%name='/home2/skurien/fractional_stats/sc1024A0001.4000.new'
%nx=1024; delx_over_eta=2.95; epsilon=3.57; teddy=1.05; % R_l=434
name='/home2/skurien/fractional_stats/sc1024A0001.4000.new'
nx=1024; delx_over_eta=2.95; epsilon=3.57; teddy=1.05; % R_l=434



ndir_use=49;
%ndir_use=49;  disp('USING ONLY 49 DIRECTIONS')

disp('1 = Scaling laws for total structure function');
disp('2 = Scaling laws and also plot D+/-');
disp('4 = 2nd and 3rd order isotropy check');
disp('5 = 2nd and 3rd order isotropy check, x,y,z directions only');
disp('6 = Scaling laws for total structure function (for paper)');
disp('7 = plot 4th order functions');
disp('8 = plot helical structure functions');
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
  [y45,y415,y43,eps,h_eps]=compisoave(name,ext,xx,ndir_use,klaws,plot_posneg,check_isotropy,1);
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
  ext=sprintf('.isostr%i%i%i',i,j,k);
  [xx,y45,y415,y43,epsl]=compisoave(name,ext,ndir_use,klaws,plot_posneg,check_isotropy);
  yysum=yysum+y45*epsl/eps/8;
  figure(10);
  plot(xx,y45);
end
end
end
  figure(10);
  plot(xx,yysum,'r');  

hold off;







