function loglog53(n,spec,time,stitle)

if (nargin==3)
  stitle=sprintf('Spectrum t=%8.4f',time);
end

x=(1:n)-1;
loglog(x,spec(1:n))
hold on
x=2:20;
y = .5*x.^(-5/3);
loglog(x,y,'r')
y = .5*x.^(-3);
loglog(x,y,'r')

axis([1,200,1e-6,1]);
title(stitle);
hold off



