%
%########################################################################
%#  plotting tacers
%########################################################################

%ts=input('time=? ');

range=0.00:1.0:1000.0 ;
%name='../src/temp';
name='/ccs/taylorm/dns/src/vxpair/vx3072_';


usefig=1;
mkpr=1;            % make ps and jpeg files


s=findstr(name,'/');
s=s(length(s));
shortname=name(s+1:length(name));

for i=range
  ts=i;
  ts = sprintf('%9.5f',10000+ts);
  ts=ts(2:10);

  fname=[name,ts,'.tracer']
  fidvor=fopen(fname,'r');
  nt=fread(fidvor,1,'float64')
  nvar=fread(fidvor,1,'float64')
  data=fread(fidvor,nt*nvar,'float64');
  if (nvar==3) 
    tracer=reshape(data,nt,nvar);
    alf=tracer(:,3);
    tracer=tracer(:,1:2);
  end

  figure(1)
  plot(tracer(:,1),tracer(:,2))
  %hold on;   plot(tracer(:,1),tracer(:,2),'r.') ;   hold off;
  axis equal
  axis([0,4,0,2]);
  title(sprintf('time = %6.2f ',i)); 

    'pause'
   % pause
end
return


