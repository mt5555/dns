%
%  2nd and 3rd order Isotropy check:  time averaged  
%
%
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

%name='/ccs/scratch/taylorm/dns/iso12/iso12_512'
%nx=512; delx_over_eta=2.74; epsilon=3.89;  teddy=1.0
%ext='.isostr001';
%times=[1:.1:7.0];


name='/ccs/scratch/taylorm/dns/sc1024A/sc1024A'
nx=2048; delx_over_eta=2.98; epsilon=3.74; teddy=1.024;
ext='.new.isostr';
times=[1:.1:1.7];


ndir_use=0;
%ndir_use=49;  disp('USING ONLY 49 DIRECTIONS')



k=0;




xx=(1:.5:(nx./2.5)) / nx;
xx_plot=(1:.5:(nx./2.5)) *delx_over_eta;   % units of r/eta

yll_ave=zeros([length(xx),15]);
ytt_ave=zeros([length(xx),15]);
ylll_ave=zeros([length(xx),15]);
yltt_ave=zeros([length(xx),15]);


times_plot=[];
for t=times
  tstr=sprintf('%10.4f',t+10000);
  fname=[name,tstr(2:10)];
  disp([fname,ext]);
  fid=fopen([fname,ext]);
  if (fid<0) ;
    disp('error openining file, skipping...');
  else
    fclose(fid);
    times_plot=[times_plot,t];
    k=k+1;
    [nx,ndelta,ndir,r_val,ke,eps_l,mu,D_ll,D_lll,D1_tt,D2_tt,D1_ltt,D2_ltt] ...
        = readisostr( [fname,ext] );
    eta_l = (mu^3 / eps_l)^.25;
    delx_over_eta_l=(1/nx)/eta_l;
    
    for dir=1:15;
      x=r_val(:,dir)/nx;                % box length
      
      yl  = D_ll(:,dir); 
      yt = .5*(D2_tt(:,dir)+D1_tt(:,dir));
      yl = spline(x,yl,xx); 
      yt = spline(x,yt,xx); 
      yll_ave(:,dir)=yll_ave(:,dir)+yl';
      ytt_ave(:,dir)=ytt_ave(:,dir)+yt';

      yl  = -D_lll(:,dir); 
      yt = -.5*(D2_ltt(:,dir)+D1_ltt(:,dir));
      yl = spline(x,yl,xx); 
      yt = spline(x,yt,xx); 
      ylll_ave(:,dir)=ylll_ave(:,dir)+yl';
      yltt_ave(:,dir)=yltt_ave(:,dir)+yt';
      
    end
  end

end
times=times_plot;
yll_ave=yll_ave/length(times);
ytt_ave=ytt_ave/length(times);
ylll_ave=ylll_ave/length(times);
yltt_ave=yltt_ave/length(times);

figure(5); clf
for i=[1]
  yll=yll_ave(:,i);
  ytt=ytt_ave(:,i)' .*xx.^(-2/3);
  ylll=ylll_ave(:,i);
  yltt=yltt_ave(:,i)' ./xx;;

  semilogx(xx_plot,ytt,'r'); hold on;
  semilogx(xx_plot,yltt,'r'); hold on;

  %
  %  compute and plot: (D_ll  + .5 r d/dr ( D_ll) )^(-2/3)
  %  
  f = yll';
  l=length(f);
  df = ( f(3:l)-f(1:l-2)) ./ (xx(3:l)-xx(1:l-2));
  f2 = f(2:l-1) + .5*xx(2:l-1).*df;
  f2 = f2 .* xx(2:l-1).^(-2/3);
  semilogx(xx_plot(2:l-1),f2,'g');

  %
  %compute and plot: [ 1/6 d/dr r D_+lll ] /r
  %
  f = ylll'.*xx/6;
  l=length(f);
  df = ( f(3:l)-f(1:l-2)) ./ (xx(3:l)-xx(1:l-2));
  df = df./xx(2:l-1);
  semilogx(xx_plot(2:l-1),df,'g');
end

axis([1 1000 0 8.0])
x=1:1000; plot(x,(4/5)*x./x,'k');
hold off;
%title('D_{lll} / r\epsilon   (4/5 law) ');
ylabel(' ','FontSize',16);
xlabel('r/\eta','FontSize',16);

print -dpsc isotime.ps

