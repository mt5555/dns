%
%########################################################################
%#  plot of DNS structure function output file
%########################################################################
%
%
clear all;

range=0:50;
%range = 1
%range=0:.5:2;

fid=fopen('../src/test.sf','r');
components=['U','V','W'];
var=['x','y','z'];

time=fread(fid,1,'float64');
while (time>=0 & time<=5000)
  disp(sprintf('plotting time=%f',time));
  nSF=fread(fid,1,'float64');           % number of SF's

  %
  % read in the U,V,W structure functions
  %
  figure(1)
  clf
  for j=1:3
  for i=1:3  % loop over x,y,z
    subplot(3,3,3*(i-1)+j)
    [n_del,delta,bin_size,n_bin,n_call,bins,pdf]=read1pdf(fid);

    if (n_call>0)
      dp=1;
      p=1;
      % structure function to the p'th power
      str=sum(pdf(:,dp).*bins.^p);

      
      mx=max(bins - bins.*(pdf(:,dp)==0));
      mn=min(bins - bins.*(pdf(:,dp)==0)); % min over non zero values      

      bar(bins,pdf(:,dp))
      ax=axis;
      %axis([-1,1,ax(3),ax(4)]);
      axis([-1.5,1.5,0,.08]);
      ylabel(['\Delta',sprintf('_{%i%s} %s',delta(dp),var(i),components(j))]);
      xlabel(sprintf('min=              %.4f  max=%.4f  ncalls=%i str(%i)=%.5f',mn,mx,n_call,p,str));
      if (i==1) title(sprintf('Time=    %.4f',time)); end;
    end
    
  end
  end

  %
  % read in the U(U**2+V**2+W**2) structure functions
  %
  figure(2)
  clf
  for j=1:3
  for i=1:3
    subplot(3,3,3*(i-1)+j)
    [n_del,delta,bin_size,n_bin,n_call,bins,pdf]=read1pdf(fid);
    
    if (n_call>0) 
      dp=1;
      p=1;
      % structure function to the p'th power
      str=sum(pdf(:,dp).*bins.^p);
      
      mx=max(bins - bins.*(pdf(:,dp)==0));
      mn=min(bins - bins.*(pdf(:,dp)==0)); % min over non zero values
      bar(bins,pdf(:,dp))
      ax=axis;
      %axis([-1,1,ax(3),ax(4)]);
      axis([-1.5,1.5,0,.05]);
      ylabel(sprintf('                  %s  DEL(%s)  p=%i',components(j),var(i),delta(dp)));
      text=['\Delta',sprintf('_{%i%s}',delta(dp),var(i))]; 
      text=[text,sprintf('              %s (',components(j)),text,'U^2 + ',text,'V^2 + ',text,'W^2)'];
      ylabel(text);
      xlabel(sprintf('min=              %.4f  max=%.4f  ncalls=%i str(%i)=%.5f',mn,mx,n_call,p,str));
      if (i==1) title(sprintf('Time=    %.4f',time)); end;
    end
    
  end
  end
  
  %
  % epsilon
  %
  % read in number of plain pdf's
  npdf=fread(fid,1,'float64');
  [n_del,delta,bin_size,n_bin,n_call,bins,pdf]=read1pdf(fid);
  figure(3)
  p=1;
  str=sum(pdf.*bins.^p);
  mx=max(bins - bins.*(pdf==0));
  mn=min(bins - bins.*(pdf==0));        % min over non zero values
  bar(bins,pdf)
  ax=axis;
  %axis([-1,1,ax(3),ax(4)]);
  axis([0,25.0,0,.01]);
  xlabel(sprintf('min=                  %.4f  max=%.4f  ncalls=%i str(%i)=%.5f',mn,mx,n_call,p,str));
  title(sprintf('Time=                  %.4f',time)); 
  
  disp('done...')
  pause  
  
  time=fread(fid,1,'float64');
end

fclose(fid);
