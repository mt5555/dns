function [f] = fun3(x,xdata)
%
%  used for fitting cubics by lsqcurvfit
%
a1 = x(1); %b1 = x(2);

f = (a1 * (xdata).^3))

