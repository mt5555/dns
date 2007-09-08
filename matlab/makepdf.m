
%
%########################################################################
%#  read in a scalar input file, plot the PDF
%
%########################################################################

name1='/tmp/balu/balu_b0040.0000.ct'
name2='/tmp/balu/balu_b0040.0000.kct'

[x,y,z]=getfield(name1);
nx=length(x);
ny=length(y);

% whole plot, but subsample data:
subsample=1;
skip=4;   x1=1; x2=nx; y1=1; y2=ny;


rangex=x1:skip:x2;
rangey=y1:skip:y2; 


[x,y,z,ct,time]=getfield(name1);
ctarray=reshape(ct,[nx*ny,1]);

qmax=max(max(max(ct)));
qmin=min(min(min(ct)));
ct = squeeze(ct(:,:,1));

if (subsample==1) 
  ct=ct(rangex,rangey); 
  x=x(rangex);
  y=y(rangey);
end

disp(sprintf('min/max ct: %f %f ',qmin,qmax));

%
% process KCT
%
[x,y,z,kct,time]=getfield(name2);
kctarray=reshape(kct,[nx*ny,1]);
qmax=max(max(max(kct)));
qmin=min(min(min(kct)));
kct = squeeze(kct(:,:,1));
if (subsample==1) 
  kct=kct(rangex,rangey); 
  x=x(rangex);
  y=y(rangey);
end
disp(sprintf('min/max kct: %f %f ',qmin,qmax));




figure(1)
subplot(2,2,1)
pcolor(x,y,ct')
shading interp
axis equal
title('ct') 

subplot(2,2,2)
hist(ctarray,40)


subplot(2,2,3)
pcolor(x,y,(kct)')
shading interp
axis equal
title('kct') 

subplot(2,2,4)
hist(kctarray./ctarray,40);






