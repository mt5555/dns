%
%########################################################################
%#  plot of DNS structure function output file
%########################################################################
%
%
clear all;


mkps=1;  % print?
dp=6;    % which seperation length index

range=0:50;
%range = 1
%range=0:.5:2;

%fid=fopen('test32.sf','r','l');
fid=fopen('iso12_256_200.sf','r','b');
%fid=fopen('../src/test32.sf','r','l');


times=[5.0 6.0 7.0 8.0 9.0 10.0];
times=[999.0];
components=['u','v','w'];
var=['x','y','z'];

n_law45=0;
time1=fread(fid,1,'float64');
while (time1>=0 & time1<=5000)
  disp(sprintf('reading time=%f',time1));
  nSF=fread(fid,1,'float64');           % number of SF's

  if (isempty(find(times==time1))) 
    plotu=0;
  else
    plotu=1;
    disp('plotting...')
  end

  if (plotu) 
    axmax=0;
    % LONGITUDINAL PDFS
    for j=1:3
      for i=1:3  
        % var(i) componenents(j)
        [n_del1,delta1,bin_size1,n_bin1,n_call1,bins1,pdf1]=read1pdf(fid);
        
        if (i==j)
          figure(j)
          subplot(4,1,1)
          bar(bins1,pdf1(:,dp))
          title(sprintf('time=%f',time1));
          ylabel(['\Delta',sprintf('_{  %i%s} %s',delta1(dp),var(i),components(j))]);
          set(gca,'YTickLabel','')
          ax=axis;
          axmax=max(axmax,ax(4));
        end
      end
    end
    
    
    % TRANSVERSE PDFS
    % components=['u','v','w'];
    % var=['x','y','z'];
    for j=1:3
      for i=1:3  
        % var(i) componenents(j)
        [n_del1,delta1,bin_size1,n_bin1,n_call1,bins1,pdf1]=read1pdf(fid);
        
        figure(i)
        subplot(4,1,j+1)
        bar(bins1,pdf1(:,dp))
        text=['\Delta',sprintf('_{%i%s}',delta1(dp),var(i))]; 
        if (i==j)
          text=[text,sprintf('%s   ||',components(i)),text,'U||^2'];
        else
          text=[text,sprintf('%s   ||',components(i)),text,components(j),'||^2'];
        end
        text=['[   ',text,'   ]^{.33}'];
        ylabel(text);
        
        set(gca,'YTickLabel','')
        ax=axis;
        axmax=max(axmax,ax(4));
        
      end
    end
    
    % make all axis the same:
    for j=1:3
      for i=1:4
        figure(j);
        subplot(4,1,i);
        axis([-4,4,0,axmax]); 
      end
    end

    if (mkps)
      figure(1);   orient tall;
      print -depsc structu.ps
      figure(2);   orient tall;
      print -depsc structv.ps
      figure(3);   orient tall;
      print -depsc structw.ps
    end
    
  else
    % longitudinal
    for j=1:3
      for i=1:3
        [n_del1,delta1,bin_size1,n_bin1,n_call1,bins1,pdf1]=read1pdf(fid);
        if (i==j) 
          m3(i)=sum(pdf1(:,dp).*bins1.^3) / ( delta1(dp)/250);
          delta_used=delta1(dp);
        end
      end
    end
    % transverse 
    for j=1:3
      for i=1:3
        [n_del1,delta1,bin_size1,n_bin1,n_call1,bins1,pdf1]=read1pdf(fid);
      end
    end
  end
  
  npdf=fread(fid,1,'float64');
  if (plotu)
    plot_epsilon(fid,time1);
  else
    [n_del1,delta1,bin_size1,n_bin1,n_call1,bins1,pdf1]=read1pdf(fid);
    epsilon=sum(pdf1(:,1).*bins1.^3);
    disp(sprintf('moment=3: %f %f %f  ave=%f  ep=%f',m3,sum(m3)/3,epsilon))
    n_law45=n_law45+1;
    law45x(n_law45)=m3(1)/epsilon;
    law45y(n_law45)=m3(2)/epsilon;
    law45z(n_law45)=m3(3)/epsilon;
    time(n_law45)=time1;
  end
  

  if (plotu)
    disp('done - press any key to continue')
    pause; 
  end;
  time1=fread(fid,1,'float64');


end
fclose(fid);

disp(sprintf('delta used = %f',delta_used));
figure(4)
plot(time,-law45x,time,-law45y,time,-law45z);
ax=axis;
axis([ax(1), ax(2), 0, 1]);

figure(5)
sum45=(law45x+law45y+law45z)/3;
plot(time,-sum45);
ax=axis;
axis([ax(1), ax(2), 0, 1]);

