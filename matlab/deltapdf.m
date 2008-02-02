%
%########################################################################
%  read and plot passive scalar pdf's
%########################################################################
%
%
clear all;


%fid=fopen('/home/mt/codes/dns/src/temp0001.0000.cpdf','r','l');
fid=fopen('../src/temp0100.0000.cpdf','r','l');

time=fread(fid,1,'float64');
npmax=fread(fid,1,'float64')         
disp(sprintf('number of pdfs = %i',npmax))
kmax = npmax/2;

figure(1); clf; 
figure(2); clf; 
for i=1:npmax
  
  [n_del,delta,bin_size,n_bin,n_call,bins,pdf]=read1pdf(fid);
  % normalize so sum(bin_size * pdf ) = 1
  % this is needed so we can plot PDFs with different binsizes on same axis:
  pdf = pdf / bin_size;

  c1=sum(bin_size * pdf.*bins);
  c2=sum(bin_size * pdf.*(bins-c1).^2);
  c4=sum(bin_size * pdf.*(bins-c1).^4);
  flatness(i)=c4/(c2^2);

  if (n_bin>350)
    % reduce bins by 50% for plotting
    count=0;
    bins2x=[];
    pdf2x=[];
    for j=2:2:n_bin
      count=count+1;
      pdf2x(count)=pdf(j)+.5*(pdf(j-1)+pdf(j+1));
      bins2x(count)=bins(j);
    end
    pdf=pdf2x; bins=bins2x;
    n_bin=count;
  end

  
  
  if (i<=npmax/2)
    k=i;
    figure(1);
    tstring = sprintf('delta filtered PDFs  k=1..%i',npmax/2);
  else
    k=i-npmax/2 ; 
    figure(2);
    tstring = sprintf('simplified delta filtered PDFs  k=1..%i',npmax/2);   
    %
  end
  disp(sprintf('k=%4i nbins=%4i  bin_size=%f',k,n_bin,bin_size));

  

  if (k==1) 
    semilogy(bins,pdf,'r*');
    title(tstring)
    hold on;
  end
  semilogy(bins,pdf);

end
figure(1); hold off;
figure(2); hold off;


figure(3); clf;
plot(1:kmax,flatness(1:kmax),'r*') ; hold on
plot(1:kmax,flatness((1:kmax) + kmax),'b*')
hold off
title('flatness delta filtered (red) and simplified (blue)')
ax=axis; axis([ax(1),ax(2),2,4]);


%orient tall
figure(1)
print('-dpng',['deltapdf','.png']); 
figure(2)
print('-dpng',['simple','.png']); 
figure(3)
print('-dpng',['flatness','.png']); 




