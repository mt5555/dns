
% points on a lighted sphere:

ndir=49;
wname=sprintf('../src/voronoi/isoave.coords%i',ndir);
w=textread(wname,'%f',1+3*2*ndir);
n=w(1);
w=w(2:length(w));
w=reshape(w,3,2*ndir);
x=w(1,:);
y=w(2,:);
z=w(3,:);


figure(1)
[xsph,ysph,zsph]=sphere(100);
surfl(xsph,ysph,zsph,'light');
colormap white
shading interp
hold on
plot3(x,y,z,'r.','MarkerSize',8);
hold off;
axis equal
print -dpsc angle.ps


