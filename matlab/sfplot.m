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

fid=fopen('test.sf','r');
components=['U','V','W'];
var=['x','y','z'];

time=fread(fid,1,'float64');
while (time>=0 & time<=5000)
  disp(sprintf('plotting time=%f',time));
  nSF=fread(fid,1,'float64');           % number of SF's
  for j=1:nSF
  figure(j)  
  for i=1:3
    subplot(3,1,i)
    % number of deltas computed
    n_del=fread(fid,1,'float64');
    delta=fread(fid,n_del,'float64');
    
    bin_size=fread(fid,1,'float64');
    n_bin=fread(fid,1,'float64');
    n_call=fread(fid,1,'float64');
    
    bins=((-n_bin:n_bin)*bin_size)';
    
    pdf=fread(fid,[(2*n_bin+1),n_del],'float64');

    dp=4;
    p=2;
    % structure function to the p'th power
    str=sum(pdf(:,dp).*bins.^p);
    
    mx=max(bins - bins.*(pdf(:,dp)==0));
    mn=min(bins - bins.*(pdf(:,dp)==0));  % min over non zero values
    bar(bins,pdf(:,dp))
    ax=axis;
    %axis([-1,1,ax(3),ax(4)]);
    axis([-2,2,0,.02]);
    ylabel(sprintf('%s(%s)-%s(%s+%i)',components(j),var(i),components(j),var(i),delta(dp)));
    xlabel(sprintf('min=%.4f  max=%.4f  ncalls=%i str(%i)=%.5f',mn,mx,n_call,p,str));
    if (i==1) title(sprintf('Time=%.4f',time)); end;
    
  end
  end
  disp('press any key...'); pause ;   disp('continuing...')
  time=fread(fid,1,'float64');
end

fclose(fid);
