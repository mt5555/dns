function  plot_epsilon(fid,time)
  %
  % epsilon
  %
  % read in number of plain pdf's
  [n_del,delta,bin_size,n_bin,n_call,bins,pdf]=read1pdf(fid);
  figure(3)
  p=3;
  str=sum(pdf.*bins.^p);
  mx=max(bins - bins.*(pdf==0));
  mn=min(bins - bins.*(pdf==0));        % min over non zero values
  bar(bins,pdf)
  ax=axis;
  axis([0,5.0,ax(3),ax(4)]);
  text='[   \mu|\nabla U|^2   ]^{1/3}';
  text=[text,sprintf('          Time=%.2f',time)];
  title(text) 
  xlabel(sprintf('[%.3f,%.3f] nc=%i <e>=%.5f',mn,mx,n_call,str));
  %title(sprintf('Time=                  %.4f',time)); 
  set(gca,'YTickLabel','')    
