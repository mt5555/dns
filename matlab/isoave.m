
name='../src/temp0000.0000';
fid=fopen([name,'.isostr'],'r','l');


ndelta=fread(fid,1,'float64');
ndir  =fread(fid,1,'float64');
nlon  =fread(fid,1,'float64');
ntran =fread(fid,1,'float64');
nnew1 =fread(fid,1,'float64');
nnew2 =fread(fid,1,'float64');

r_val=fread(fid,[ndelta,ndir],'float64');
if (nlon==2) 
   D_ll=fread(fid,[ndelta,ndir],'float64');
   D_lll=fread(fid,[ndelta,ndir],'float64');
end
if (ntran==2) 
   D1_tt=fread(fid,[ndelta,ndir],'float64');
   D2_tt=fread(fid,[ndelta,ndir],'float64');
   D1_ltt=fread(fid,[ndelta,ndir],'float64');
   D2_ltt=fread(fid,[ndelta,ndir],'float64');
end


figure(1)
for i=1:ndir
  semilogx( r_val(:,i),D_ll(:,i) )
  hold on;
    semilogx( r_val(:,i),D_ll(:,i),'o' )
end
hold off;

figure(2)
for i=1:ndir
  semilogx( r_val(:,i),D_lll(:,i) )
  hold on;
    semilogx( r_val(:,i),D_lll(:,i),'o' )
end
hold off;


