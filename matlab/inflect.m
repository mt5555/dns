function[indy] = inflect(x,y,color)
     c = ['b','k','r','c','g','y'];
     len = length(x)-1;
     d2 = diff(diff(log10(y))./diff(log10(x)))./diff(log10(x(1:len)));
     semilogx(x(1:len-1),d2,c(color))
     indx = find(x > 20 & x < 800);
%     max(d2(indx));
%     d2(indx);
% index find for integer str fns 1024^3
%     indy = find(d2(indx) <= max(d2(indx)) & d2(indx) >= max(d2(indx))-0.4);

% index find for for fractional str fns 1024^3

     indy  = find(d2(indx) <= max(d2(indx)) & d2(indx) >= max(d2(indx))-0.05);

%shift back to proper x-axis index;
     indy = indy + indx(1); 

