% sphereplot
% points on a lighted sphere:

ndir=73;
wname=sprintf('../src/voronoi/isoave.coords%i',ndir);
w=textread(wname,'%f',1+3*2*ndir);
n=w(1);
w=w(2:length(w));
w=reshape(w,3,2*ndir);
x=w(1,:);
y=w(2,:);
z=w(3,:);


figure(1)
[xsph,ysph,zsph]=sphere(200);
surfl(xsph,ysph,zsph,'light');
colormap white
shading interp
hold on
scale=1.01;
%plot3(scale*x,scale*y,scale*z,'k.','MarkerSize',12);
plot3(scale*x,scale*y,scale*z,'r.','MarkerSize',12);
hold off;
axis equal
text(0,-1.4,-1.1,'X')
text(-1.5,0,-1.1,'Y')
text(-1.3,1.2,.1,'Z')

print -deps  angle.ps
print -dtiff angle.tiff


