function loglog53(n,spec,time,stitle)

if (nargin==3)
  stitle=sprintf('Spectrum t=%8.4f',time);
end

s=size(spec);
s=s(2);
ccol=[ 'b','g','r' ];  


x=(1:n)-1;

scale = 1;
%scale=.1*x.^(5/3);

for i=1:s
   loglog(x, scale.*spec(1:n,i)',ccol(i))
   hold on
end

x=2:300;
scale = 1;
%scale=x.^(5/3);

y = .5*x.^(-5/3);
loglog(x,scale.*y,'k')

%y = .5*x.^(-3);
%loglog(x,scale.*y,'k')

axis([1,1e4,1e-7,1]);
title(stitle);
hold off



