%########################################################################
%
% read in and plot velocity increment PDFs 
%
% (see also more complex and older sfplot.m)
%
%########################################################################
%
%
clear all;


dp=2;    % which seperation length index to plot

fid=fopen('/home/mataylo/codes/dns/src/temp0002.0000.sf','r','l');

components=['u','v','w'];
var=['x','y','z'];

  time1=fread(fid,1,'float64');
  disp(sprintf('reading time=%f',time1));
  nSF=fread(fid,1,'float64');           % number of PDFs in file

  % LONGITUDINAL PDFS  
  for j=1:3
    for i=1:3  
      % var(i) componenents(j)
      [n_del1,delta1,bin_size1,n_bin1,n_call1,bins1,pdf1]=read1pdf(fid);
      % ignore the i<>j PDFs
      if (i==j)
        subplot(3,1,i)
        semilogy(bins1(:,dp),pdf1(:,dp))
        title(sprintf('time=%f',time1));
        ylabel(['\Delta',sprintf('_{  %i%s} %s',delta1(dp),var(i),components(j))]);
        set(gca,'YTickLabel','')
      end
    end
  end
  
  % sometimes higher order PDFs are also in the input file,
  % skip them.  we've read 3 PDFs so far, so skip the rest:
  for i=[4:nSF]
    for j=1:3
      % these are the 3rd and 4th order SF 
      [n_del1,delta1,bin_size1,n_bin1,n_call1,bins1,pdf1]=read1pdf(fid);
    end
  end

% final PDF for each snapshot is for the scalar epsilon 
npdf=fread(fid,1,'float64');
% plot_epsilon(fid,time1);  % read and plot epsilon pdf


fclose(fid);

