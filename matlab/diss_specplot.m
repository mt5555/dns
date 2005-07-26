%read and compute mean D_E
%D_E = 2*nu*Sum_0^k (p^2 E(p))*(2*pi)^2


namedir = '/nh/nest/u/skurien/projects/helicity_data/helical_forced/hel512_hpi2/diag/';
name = 'skhel512_hpi2';
mu = 1e-4;

pname = [strrep(name,'_','-'),'.spec'];
fid=fopen([namedir,name,'.spec'],'r','l');
disp([namedir,name,'.spec']);

movie=0;

time=fread(fid,1,'float64');
time
spec_r = [];
spec_ave = spec_r;

j=0;l=0;
while (time>=.0 & time<=9999.3)
     n_r=fread(fid,1,'float64');
     spec_r=fread(fid,n_r,'float64');
     %if (j==0)
     if (time < 1);
     spec_ave =  0*spec_r;  %initialize spec_ave
     end
     if (time >= 1);        %only average times greater than 1 eddy turnover
     j=j+1;                 %count js for averaging only after t==1.
     spec_ave = spec_ave + spec_r;
     end
     j
     
n_x=fread(fid,1,'float64');
spec_ux=spec_scale*fread(fid,n_x,'float64');
spec_vx=spec_scale*fread(fid,n_x,'float64');  
spec_wx=spec_scale*fread(fid,n_x,'float64');  

n_y=fread(fid,1,'float64');
spec_uy=spec_scale*fread(fid,n_y,'float64');
spec_vy=spec_scale*fread(fid,n_y,'float64');
spec_wy=spec_scale*fread(fid,n_y,'float64');

n_z=fread(fid,1,'float64');
spec_uz=spec_scale*fread(fid,n_z,'float64');
spec_vz=spec_scale*fread(fid,n_z,'float64');
spec_wz=spec_scale*fread(fid,n_z,'float64');  

   time=fread(fid,1,'float64');
   time
end
%compute total Energy for snapshot
E = sum(spec_r);
disp(sprintf('Total Energy  = %d',E));

k = [0:n_r-1];
%compute dissipation rate of energy
k=k'
%epsilon = mu*sum(spec_r .* (k*2*pi).^2);
disp(sprintf('energy dissipation rate = %d',epsilon));

end

fclose(fid);  

%end

% mean energy spectrum
spec_ave = spec_ave/j;
spec_l = spec_l/l;
spec_t1 = spec_t1/l;
spec_t2 = spec_t2/l;
spec_t3 = spec_t3/l;
spec_t = spec_t/l;

%compute mean dissipation rate of energy
epsilon = 2*mu*sum(spec_ave .* (k*2*pi).^2);
disp(sprintf('mean energy dissipation rate = %d',epsilon));


for i= 0:n_r-1
D_E(i+1) = 2*mu*sum(spec_ave(1:i+1).* (k(1:i+1)*2*pi).^2);
end

figure(5);clf;
loglog(k,abs(D_E - epsilon)/epsilon,'k--'); hold on;
title('Mean energy dissipation spectrum');



