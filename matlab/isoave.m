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


%name='/ccs/scratch/taylorm/check256_0000.8000'
%name='/ccs/scratch/taylorm/dns/iso12_5120002.7000'
name='/ccs/scratch/taylorm/dns/iso12_5120003.0000'
%ext='.isostr001';



ndir_use=0;
%ndir_use=49;  disp('USING ONLY 49 DIRECTIONS')




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
   ndir_use=3;
   klaws=0;
   check_isotropy=1;
else
   klaws=1;
   check_isotropy=0;
   if (in==2) 
     plot_posneg=1;
   end
end



[xx,y45,y415,y43,eps]=compisoave(name,ext,ndir_use,klaws,plot_posneg,check_isotropy);

figure(10); hold off;
semilogx(xx,y45,'k'); hold on;
axis([1,1000,0,1.5]);



return



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







