function loglog53(n,spec,stitle,CK,option)


ccol=[ 'b','g','r','c','m','y', 'b','g','r','c','m','y' ];  

if (nargin==4) option=0; end;

% default, if nothing specified
if (option==0) 
  x_ck=1:500;
  % x_ck=-1;  % disable scale line
   scalep = 0;
   %scalep = 5/3;
   
   plaw=-5/3;
   %plaw=-3 

   ax=[1,1e4,1e-6,1.0];
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
   scalep = 0;
   plaw=-5/3;
   ax=[1,1e4,1e-6,1.0];
end


x=(1:n)-1;
scale=x.^scalep;
s=size(spec);
s=s(2);
for i=1:s
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



