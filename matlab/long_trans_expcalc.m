% plot longitudinal and transverse structure functions and compute their
% scaling exponents

function[] = long_trans_expcalc(Dl,Dt,ndelta,ndir,r_val,nx,delx_over_eta,xx)

%for i = 1:2:9
% for j = 1:ndir
%figure(10)
%     loglog(r_val(:,j), abs(Dl(:,j,i))); hold on;
%figure(11)
%     loglog(r_val(:,j),abs(Dt(:,j,i)));hold on;
%end
%end

%perform angle-average

if (ndir==3)
  w=ones([1,3])/3;
else
  equalw=0;
  if (equalw) 
    % put this in to use equally weighted:
    w=ones([1,ndir])/ndir;
  else
    % get the weights:
    wname=sprintf('../src/voronoi/isoave.weights%i',ndir);
    w=textread(wname,'%f');
    % take every other weight
    w=2*w(1:2:length(w));
  end
end
if (abs(1-sum(w))>1e-7) 
  disp('error: weights do not sum to 1')
  return;
end

yyavel=0*xx;
yyavet=0*xx;
for i = 1:2:5
for dir=1:ndir
  x = r_val(:,dir);                       % units of box length
  x_plot=x*nx*delx_over_eta;  % units of r/eta
  xx_plot = xx*nx*delx_over_eta;        % units of r/eta
  yl = Dl(:,dir,i);
  yt = Dt(:,dir,i);
  yyavel=yyavel+w(dir)*spline(x,yl,xx);
  yyavet = yyavet+w(dir)*spline(x,yt,xx);
end

  figure(12)
  loglog(xx_plot,yyavel,'r'); hold on;
  loglog(xx_plot,yyavet,'k'); hold on;
xlabel('r/\eta');
    title('Structure functions - longitudinal (red), transverse (black)');

%LONGITUDINAL AND TRANSVERSE LOCAL SLOPES
figure(13)
 derl = localslp(xx_plot,yyavel,3);hold on;
 dert = localslp(xx_plot,yyavet,2);hold on
 title('Local slopes - longitudinal (red), transverse (black)');

% LONGITUDINAL 
 
 %extra stuff for fractional statistics paper
 if(0)
 figure(15)
 ind = inflect(xx_plot,yyavel,3);hold on;
 title('Inflection points - longitudinal');
 p=i+1;
 sprintf('p=%d;',p)
 sprintf('Exps = %f; ',mean(derl(ind)))
 sprintf('Error on exps = %f',std(derl(ind)))
     figure(17)
     k41 = p/3;
     plot(p,mean(derl(ind)),'o',p,k41,'x');hold on;
     title('Comparison of exps with K41 - longitudinal');
plot(k41,mean(derl(ind)),'ro'); hold on;

figure(19)
     plot(p,(mean(derl(ind))-k41)/k41,'ko');hold on;
     title('Relative scaling exponents - longitudinal');

     figure(21)
     plot(p,mean(derl(ind)),'o');hold on;
     title('comparison of long and trans exps');
 end

% TRANSVERSE
%extra stuff for fractional statistics paper
if(0)
figure(16)
ind = inflect(xx_plot,yyavet,3);hold on;
title('Inflection points - transverse');

sprintf('p=%d;',p)
sprintf('Exps = %f; ',mean(dert(ind)))
sprintf('Error on exps = %f',std(dert(ind)))
     figure(18)
     k41 = p/3;
     plot(p,mean(dert(ind)),'o',p,k41,'x');hold on;
     title('Comparison of exps with K41 - transverse');
%plot(k41,mean(dert(ind)),'ro'); hold on;

figure(20)
     plot(p,(mean(dert(ind))-k41)/k41,'ko');hold on;
     title('Relative scaling exponents -  transverse');

     figure(21)
     plot(p,mean(dert(ind)),'x');hold on;
     title('comparison of long and trans exps');
end

end

xlabel('r/\eta');
%ax=axis;  axis([1,xmax,ax(3),ax(4)]);
%ax=axis;  axis([1,xmax,0,.3]);
%hold off;
%if (plot_points==1) 
%print('-dpsc',[bname,'_415.ps']);
%print -djpeg 415.jpg
