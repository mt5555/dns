function[d1] = localslp(x,y,color)
     c = ['b','k','r','c','g','y'];
     len = length(x)-1;
     d1 = diff(log10(y))./diff(log10(x));
%    semilogx(x(1:len),d1,c(color))
