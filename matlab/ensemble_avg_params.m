function [avg_eps, avg_heps, avg_delx_over_eta] = ensemble_avg_params(name,ext,times)


avg_eps = 0;
avg_heps = 0;
avg_delx_over_eta = 0;

for t=times;
  tstr=sprintf('%10.4f',t+10000);
  fname=[name,tstr(2:10)];
  disp([fname,ext]);
  fid=fopen([fname,ext]);
  if (fid<0) ;
    disp('error openining file, skipping...');
    else
    fclose(fid);
  end
[nx,ndelta,ndir,r_val,ke,eps_l,mu,tmp,tmp,tmp,tmp,tmp,tmp,tmp,...
    tmp,tmp,tmp,tmp,tmp,H_ltt,H_tt,tmp,tmp,tmp,h_eps_l] ...
        = readisostr( [fname,ext] );

eta = (mu^3 /eps_l)^.25;
delx_over_eta=(1/nx)/eta;

avg_delx_over_eta = avg_delx_over_eta + delx_over_eta/30;
avg_eps = avg_eps + eps_l/30;
avg_heps = avg_heps + h_eps_l/30;

end
