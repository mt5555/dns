%read and compute mean k*E(k). Use to calculate relative helicity in conjunction%with hspecplot

%namedir = '/scratch2/taylorm/sc1024A/';
%name = 'sc1024A';
%mu=.35e-4;

namedir = '/nh/nest/u/skurien/projects/helicity_data/helical_forced/hel512_hpi2/diag/';
name = 'skhel512_hpi2';
mu = 1e-4;

%namedir = '/nh/nest/u/skurien/projects/helicity_data/helical_forced/hel512_hpi2/diag/';
%name = 'sc1024A';
%mu = 3.5e-5;

% plot all the spectrum: 
movie=0;

% plot the relative helicity 
relh = 1;


pname = [strrep(name,'_','-'),'.kspec'];
fid=fopen([namedir,name,'.kspec'],'r','l');


kspec = [];
kspec_ave = kspec;

[time,count]=fread(fid,1,'float64');
j = 0;
while (time >=.0 & time <= 9999.3)
if (count==0) 
   disp('error reading kspec file')
end
time
n_r=fread(fid,1,'float64');
kspec = fread(fid,n_r,'float64');

if (time < 1)
     kspec_ave = 0*(kspec); % initialize kspec_ave
end

if (time >= 1)
j=j+1
if (length(kspec_ave)==0) 
  kspec_ave = kspec;
else
  kspec_ave = kspec_ave + kspec;
end
end

[time,count]=fread(fid,1,'float64');
end

kspec_ave = kspec_ave/(j);

%plot relative helicity hspec_ave/2/kspec_ave
if (relh ==1)
     figure(13);clf
semilogx(k,abs(hspec_ave)./(kspec_ave),'k'); % kspec_ave already has factor 2

end
