%
%########################################################################
%  read and plot passive scalar pdf's
%########################################################################
%
%
clear all;


fid=fopen('/scratch2/taylorm/tmix256C/tmix256C0001.0001.spdf','r','l');

time1=fread(fid,1,'float64');
npmax=fread(fid,1,'float64');         


np=1;
for p=1:npmax
   [n_del,delta,bin_size,n_bin,n_call,bins,pdf]=read1pdf(fid);
   if (p==np) break; end;
end



figure(1); clf; subplot(1,1,1)
s2=sum(pdf.*bins.^2);
mx=max(bins - bins.*(pdf==0));
mn=min(bins - bins.*(pdf==0));        % min over non zero values
semilogy(bins,pdf)
%  ax=axis;
%  axis([0,5.0,ax(3),ax(4)]);
xlabel(sprintf('[%.3f,%.3f] <s^2>=%.5f',mn,mx,s2));
title(sprintf('passive scalar  t=%.4f',time)); 
set(gca,'YTickLabel','')    
