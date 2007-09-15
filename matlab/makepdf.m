
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
subsample=0;
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

disp(sprintf('min/max cos(theta): %f %f ',qmin,qmax));

%
% process KCT
%
[x,y,z,kct,time]=getfield(name2);
kappa=kct./ct;
kctarray=reshape(kct,[nx*ny,1]);
kappa_array=reshape(kappa,[nx*ny,1]);
qmax=max(max(max(kct)));
qmin=min(min(min(kct)));
kct = squeeze(kct(:,:,1));
if (subsample==1) 
  kct=kct(rangex,rangey); 
  x=x(rangex);
  y=y(rangey);
end
disp(sprintf('min/max k cos(theta): %f %f ',qmin,qmax));


pr=1:64;

figure(1)
subplot(2,2,1)
pcolor(x(pr),y(pr),ct(pr,pr)')
colorbar
shading flat
axis equal
title('cos(\theta)  (lower left corner)') 

subplot(2,2,2)
[h,xh]=hist(acos(ctarray),100);
h=h/sum(h);
xh=xh/pi;
plot(xh,h);
axis([0 1 0.008 .012]);
xticks(['    0';' pi/4';' pi/2';'3pi/4';'   pi'],' ');
grid
title('PDF of \theta ')


subplot(2,2,3)
pcolor(x(pr),y(pr),(kct(pr,pr))')
caxis([-.001,.001]);
colorbar
shading flat
axis equal
title('\kappa cos(\theta)  (lower left corner)') 

% $$$ kctp = (kctarray>0).*kctarray;
% $$$ kctp = (kctp==0)*1e-12 + kctp;
% $$$ min(kctp)
% $$$ kctn = (kctarray<0).*kctarray;
% $$$ kctn = -(kctn==0)*1e-12 + kctn;
% $$$ max(kctn)


subplot(2,2,4)
[h,xh]=hist(kctarray,200);
h=h/sum(h);
[h2,xh2]=hist(kappa_array,200);
h2=h2/sum(h2);

semilogy(xh,h,xh2,h2)
title('PDF      \kappa cos(\theta) (blue)       \kappa (green)')
axis([-.01 .01 1e-6 10 ]);
grid
print -r120 -dpng pdfkappa.png 




