function loglog53(n,spec,stitle,CK)




s=size(spec);
s=s(2);
ccol=[ 'b','g','r','c','m','y', 'b','g','r','c','m','y' ];  

%scalep = 0;
scalep = 5/3;

x=10:500;
scale=x.^scalep;



y = CK* scale .* x.^(-5/3);
%y = CK* scale .* x.^(-3);
%semilogx(x,y,'k'); hold on
loglog(x,y,'k'); hold on


x=(1:n)-1;
scale=x.^scalep;
for i=1:s
   plot(x, scale.*spec(1:n,i)',ccol(i))
end


if (scalep==0) 
%   axis([1,1e4,1e-6,.1]);
   axis([1,200,1e-12,10.0]);
else
   axis([1,1e4,1e-5,1e1]);
%   axis([1,1e4,0,1.0]);
end
title(stitle);
hold off



