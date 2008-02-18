function loglog53(n,spec,stitle,CK,option)
%
%
%  plot spectrum  EP(k)  = scale1 * E(k) * k^(scalep)
%  
%  if length(x_ck)>1, plot a slope line:    CK* x_ck.^plaw;
%
%  option is used to specify sets of (scalep,scale1,x_ck)
%

ccol=[ 'b','g','r','c','m','y', 'b','g','r','c','m','y' ];  
x=(1:n)-1;
scale1=1;

if (nargin==4) option=0; end;

% default, if nothing specified
if (option==0) 
  x_ck=1:500;
  % x_ck=-1;  % disable scale line
   plaw=-5/3;

   scalep = 0;
   ax=[1,1e3,1e-6,1.0];

   %scalep = 5/3;   % scale out for a 5/3 slope
   %ax=[1,1e3,1e-2,10.0];
                      
   %plaw=-3 
   %ax=[1,1e2,1e-6,1.0];
   

end

% decay2048 KE:
if (option==1) 
   x_ck=1:500;
   scalep = 0;
   plaw=-5/3;
   ax=[1,1e4,1e-7,1.0];
end
% decay2048 enstropy:
if (option==2) 
   x_ck=1:500;
   scalep = 0;
   plaw=1/3;
   ax=[1,1e4,1e-5,100.0];
end

% used by passive scalars
if (option==3) 
   x_ck=1:500;
   plaw=-5/3;
   scalep = 0;
   ax=[1,1e4,1e-12,1e-4];
   %scalep = 4;
   %ax=[1,1e4,1e-9,1000.0];
end

% used by passive scalars, 2D
if (option==4) 
   x_ck=1:1000;
   plaw=-3;
   scalep = 0;
   ax=[1,1e4,1e-6,1.0];
   ax=[1,1e3,1e-12,1e-2];
end
% decay2048 KE:
if (option==5) 
   x_ck=1;      % dont plot 5/3 line
   scalep = 0;
   plaw=-5/3;
   ax=[1,1e4,1e-7,1.0];
   x=1.25*x;    % rescale x axis so peak match KCM plot
   scale1=3;    % rescale E to match KCM plot 
end
% used by aspect ratio 2D spectrum
if (option==6) 
   x_ck=10:100;     
   scalep = 0;
   plaw=-3;
   ax=[1,1e2,1e-10,1e-1];
end
if (option==7) 
   x_ck=10:100;     
   scalep = 0;
   plaw=2;
   ax=[1,1e3,1e-6,1];
end



scale=scale1*x.^scalep;
s=size(spec);
s=s(2);
for i=1:s
%   loglog(x, scale.*spec(1:n,i)',ccol(i),'LineWidth',2.0); hold on;
   loglog(x, scale.*spec(1:n,i)',ccol(i)); hold on;
end

if (length(x_ck)>1)
  scale=x_ck.^scalep;
  y = CK* scale .* x_ck.^plaw;
  loglog(x_ck,y,'k'); 
end

axis(ax);

title(stitle);
hold off



