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
name='/ccs/scratch/taylorm/dns/iso12_512'
nx=512; delx_over_eta=2.74; epsilon=3.89;
%ext='.isostr001';



ndir_use=0;
%ndir_use=49;  disp('USING ONLY 49 DIRECTIONS')



k=0
times=[0:.1:3.7];
for t=times
  tstr=sprintf('%10.4f',t+10000);
  fname=[name,tstr(2:10)];
  k=k+1;

  klaws=1;  % compute 4/5 laws
  plot_posneg=0;
  check_isotropy=0;
  [xx,y45,y415,y43,eps]=compisoave(fname,ext,ndir_use,klaws,plot_posneg,check_isotropy);

  mx45_localeps(k)=max(y45);
  mx45(k)=max(y45)*eps/epsilon;
  
end
end
end


figure(4); hold off;
plot(times,mx45,times,mx45_localeps);








