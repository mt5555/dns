%
%########################################################################
%#  plotting tacers
%########################################################################

%ts=input('time=? ');

range=0:.1:50.0 ;
%name='../src/temp';
%name='/ccs/taylorm/dns/src/vxpair/vx3072_';
%name='/data/vxpair/vx2048c';
name='/ccs/taylorm/dns/src/vxpair/vx4096d';


usefig=1;
mkpr=1;            % make ps and jpeg files

ccol=[ 'b','g','r','c','m','y', 'b','g','r','c','m','y' ];  


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
    
    figure(1); clf;
    plot(tracer(1:nt_ins,1),tracer(1:nt_ins,2))
    hold on;   
    plot(tracer(nt_ins+1:nt,1),tracer(nt_ins+1:nt,2),'r.');
    %for k=1:8
    %  plot(vxline_x(k,:),vxline_y(k,:),ccol(k));
    %end 
    hold off;
    axis equal
    axis([0,4,0,2]);
    title(sprintf('time = %6.2f ',i)); 
%    'pause' ;     pause
  end
end

figure(2)

clf;
  hold on;   
  for k=1:8
    plot(vxline_x(k,:),vxline_y(k,:),ccol(k));
    plot(vxline_x(k,1),vxline_y(k,1),'k.');
  end 

  hold off;
  axis image
  title(sprintf('time = %6.2f ',time)); 

return


