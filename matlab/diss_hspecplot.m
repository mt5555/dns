%read and compute mean D_H
%D_H = 2*nu*Sum_0^k (p^2 H(p))*(2*pi)^2


namedir = '/nh/nest/u/skurien/projects/helicity_data/helical_forced/hel512_hpi2/diag/';
name = 'skhel512_hpi2';
mu = 1e-4;

pname = [strrep(name,'_','-'),'.hspec'];
fid=fopen([namedir,name,'.hspec'],'r','l');
disp([namedir,name,'.hspec']);

movie=0;


hspec_n = [];
hspec_p = [];
hspec_ave = hspec_n + hspec_p;

[time,count]=fread(fid,1,'float64');
time

j=0;
while (time>=.0 & time<=9999.3)
     n_r=fread(fid,1,'float64');
     hspec_n=fread(fid,n_r,'float64');
     hspec_p=fread(fid,n_r,'float64');
     %if (j==0)
     if (time < 1);
     hspec_ave =  0*(hspec_n + hspec_p);  %initialize hspec_ave
     end
     if (time >= 1);        %only average times greater than 1 eddy turnover
     j=j+1;                 %count js for averaging only after t==1.
     if (length(hspec_ave)==0) 
     hspec_ave = hspec_n + hspec_p;
     else
     hspec_ave = hspec_ave + hspec_n + hspec_p;
end
     end
     j
% compute total helicity for snapshot -- in Mark's units
% multiply by 2*pi^2 to get Takeshi's units
H = sum(hspec_p+hspec_n);
disp(sprintf('Total helicity = %d',H))

% compute dissipation rate of helicity -- in Mark's units 
% multiply dissipation rate by 2*pi to get Takeshi's units

h = -2*mu*sum((hspec_p+hspec_n).*(k*2*pi).^2);
disp(sprintf('Helicity dissipation rate = %d',h));

[time,count]=fread(fid,1,'float64');
     

end

%compute total Energy for snapshot
E = sum(spec_r);
disp(sprintf('Total Energy  = %d',E));

k = [0:n_r-1];
%compute dissipation rate of helicity
k=k';
h = mu*sum((hspec_p+hspec_n) .* (k*2*pi).^2);
disp(sprintf('helicity dissipation rate = %d',epsilon));

end

fclose(fid);  
hspec_ave = hspec_ave/(j);

%end


%compute mean dissipation rate of helicity
h = 2*mu*sum(hspec_ave .* (k*2*pi).^2);
disp(sprintf('mean helicity dissipation rate = %d',h));


for i= 0:n_r-1
D_H(i+1) = 2*mu*sum(hspec_ave(1:i+1).* (k(1:i+1)*2*pi).^2);
end

figure(5);
loglog(k,(h- D_H)/h,'r'); hold on;
title('Mean helicity dissipation spectrum');



