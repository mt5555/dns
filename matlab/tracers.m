%########################################################################
%#  plotting tacers
%########################################################################
usefig=1;
mkpr=0;            % make ps and jpeg files
mkascii=0;         % output ascii point data
plot_noni=1;       % plot the non-insert points also
axis_size=[0,4,0,2];
ccol=[ 'b','g','r','c','m','y', 'b','g','r','c','m','y' ];  

%ts=input('time=? ');

range=0:5:200.0 ;
range=0:.1:4.0;
%name='../src/temp';
%name='/ccs/taylorm/dns/src/vxpair/vx3072_';
%name='/data/vxpair/vx2048c';
%name='/ccs/taylorm/dns/src/vxpair/vx4500a';
%name='/ccs/taylorm/dns/src/vxpair/vx4096d';
%name='/scratch2/taylorm/vx12288b/vx12288b';
%name='/home/scratch/vxpair/vx12288/vx12288b';



%name='/ccs/scratch/taylorm/kras/vx2560a/vx2560a'; 
%name='/ccs/scratch/taylorm/kras/vx2560b/vx2560b'; 
%name='/ccs/scratch/taylorm/kras/vx2560c/vx2560c'; 
%name='/ccs/scratch/taylorm/kras/vx5120a/vx5120a';
%name='/home/scratch/kras/vx1280a/vx1280a'; 
%plot_noni=0;   mkpr=1;  axis_size=[0,2.5,0,2];

name='/home/scratch/kras/vx7680a/vx7680a'
plot_noni=0;   mkpr=1;  axis_size=[.4,1.8,0,1.1];





s=findstr(name,'/');
s=s(length(s));
shortname=name(s+1:length(name));
count=0;

for i=range
  ts=i;
  ts = sprintf('%9.5f',10000+ts);
  ts=ts(2:10);

  fname=[name,ts,'.tracer']
  fidvor=fopen(fname,'r');
  if (fidvor>=0) 
    nt=fread(fidvor,1,'float64');  % total number of points including insert
    nt_ins=fread(fidvor,1,'float64');      % total number of non-insert points
    nvar=fread(fidvor,1,'float64')
    data=fread(fidvor,nt*nvar,'float64');
    count=count+1;
    time=i;
    if (nvar==3) 
      tracer=reshape(data,nt,nvar);
      alf=tracer(:,3);
      tracer=tracer(:,1:2);
    end
    
    if (i==range(1))
      vxline_x=tracer(nt_ins+1:nt,1);
      vxline_y=tracer(nt_ins+1:nt,2);
    else
      vxline_x=[vxline_x , tracer(nt_ins+1:nt,1)];
      vxline_y=[vxline_y , tracer(nt_ins+1:nt,2)];
    end

    if (1)     
      figure(1); clf;
      plot(tracer(1:nt_ins,1),tracer(1:nt_ins,2))
      if (plot_noni)
        hold on;   
        plot(tracer(nt_ins+1:nt,1),tracer(nt_ins+1:nt,2),'r.');
        hold off;
      end
      axis equal
      axis(axis_size);
      title(sprintf('time = %6.2f ',i)); 
      if (mkpr)
         pname=[name,'.tracer.ps'];
         print('-depsc',pname);
      end
      'pause' ;     pause
    end
    
    % output to ASCI file
    if (mkascii) 
      fidout=fopen('ascii.out','w'); 
      for i=1:nt_ins
        fprintf(fidout,'%20.12f %20.12f\n',tracer(i,1),tracer(i,2));
      end
      fclose(fidout);
    end
  end
end

if (plot_noni)
  figure(2);
  clf;
  hold on;   
  for k=1:8
    plot(vxline_x(k,:),vxline_y(k,:),ccol(k));
    plot(vxline_x(k,1),vxline_y(k,1),'k.');
  end 
  
  hold off;
  axis image
  title(sprintf('time = %6.2f ',time)); 
end
return


