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

fid=fopen('test64.sf','r');


times=[1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0 10.0];
mkps=1;
dp=5;

time=fread(fid,1,'float64');
while (time>=0 & time<=5000)
  disp(sprintf('reading time=%f',time));
  nSF=fread(fid,1,'float64');           % number of SF's

  if (isempty(find(times==time))) 
    plotu=0;
    plottran=0;
    plote=0;
  else
    plotu=1;
    plottran=1;
    plote=1;
    disp('plotting...')
  end
  
  if (plotu) 
    figure(1)
    plot_u(fid,time,dp);
    if (mkps) 
      orient landscape
      print -dpsc sf1.ps
    end
  else
    % skip data
    for i=1:9; tmp=read1pdf(fid);  end;
  end

  if (plottran)
    figure(2)
    plot_tran(fid,time,dp);
    if (mkps) 
      orient landscape
      print -dpsc sf2.ps
    end
  else
    % skip data
    for i=1:9, tmp=read1pdf(fid); end;
  end

  npdf=fread(fid,1,'float64');
  if (plote)
    plot_epsilon(fid,time);
  else
    tmp=read1pdf(fid);
  end
  

  if (plotu | plottran | plote) 
    disp('done - press any key to continue')
    pause; 
  end;
  time=fread(fid,1,'float64');
end

fclose(fid);
