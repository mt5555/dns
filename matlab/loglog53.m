function loglog53(n,spec,time)

x=(1:n)-1;
loglog(x,spec)
hold on
x=2:20;
y = .5*x.^(-5/3);
loglog(x,y,'r')
axis([1,100,1e-5,1]);
title(sprintf('Spectrum t=              %8.4f',time));
hold off



