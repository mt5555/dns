%
%########################################################################
%  read and plot scalar pdf's ( *.cpdf files )
%########################################################################
%  flatness(k=1)  2.7050   2.59
%                 2.7500   2.68
clear all;


%fid=fopen('/home/mataylo/tmp/hyp5120002.7500.cpdf','r','l');
fid=fopen('/home/mataylo/tmp/hyp5120002.7050.cpdf','r','l');
%fid=fopen('./sc1024A0002.5000.cpdf','r','l');

time=fread(fid,1,'float64');
npmax=fread(fid,1,'float64')         
disp(sprintf('number of pdfs = %i',npmax))
kmax = npmax;

figure(1); clf; 
for k=1:npmax
  
  [n_del,delta,bin_size,n_bin,n_call,bins,pdf]=read1pdf(fid);
  % normalize so sum(bin_size * pdf ) = 1
  % this is needed so we can plot PDFs with different binsizes on same axis:
  pdf = pdf / bin_size;
  c1=sum(bin_size * pdf.*bins);
  c2=sum(bin_size * pdf.*(bins-c1).^2);
  c4=sum(bin_size * pdf.*(bins-c1).^4);
  flatness(k)=c4/(c2^2);
  

  if (n_bin>3500)
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

  
  
  figure(1);
  tstring = sprintf('delta filtered PDFs  k=1..%i',npmax);

  if (k==1) 
    semilogy(bins,pdf,'r*');
    title(tstring)
    hold on;
  end
  if (k<16) 
     disp(sprintf('k=%4i nbins=%4i  bin_size=%e flatness=%e',k,n_bin,bin_size,flatness(k)));
     plot(bins,pdf);
  end
end
figure(1); hold off;


figure(2); clf;
plot(1:kmax,flatness(1:kmax),'r*') ; hold on
hold off
title('flatness delta filtered')
ax=axis; axis([ax(1),ax(2),1,4]);


%orient tall
figure(1)
print('-dpng',['deltapdf','.png']); 
figure(2)
print('-dpng',['flatness','.png']); 




