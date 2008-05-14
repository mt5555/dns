%
%########################################################################
%  read in a famil of ".cpdf" files all with the same binsize
%  and merge them into 1 pdf
%########################################################################
clear all;


% base name. time will be appended for each file
name = '/home/mataylo/tmp/hyp512';

% merge these times:
times = 2.7050:.0050:2.75;

count = 0;
maxbin = 0;

% go through all PDFs, finding the maximum number of bins
% and verify binsizes for each k are always the same
disp('scanning all PDFs...')
for t = times
  tstr=sprintf('%10.4f',t+10000);
  fname=[name,tstr(2:10),'.cpdf'];

  fid=fopen(fname,'r','l');
  if (fid>0) 
    disp(fname);
    count = count + 1;

    time=fread(fid,1,'float64');
    npmax=fread(fid,1,'float64');         
    for k=1:npmax
      [n_del,delta,bin_size,n_bin,n_call,bins,pdf]=read1pdf(fid);
      maxbin = max(maxbin,n_bin);

      if (count==1) 
        bin_size_list(k)=bin_size;
      else
        if (bin_size_list(k) ~= bin_size) then
          disp(sprintf('Error: bin sizes do not match: k=%i  bin sizes %e %e ',k,bin_size_list(k),bin_size))
        end
      end
    end
    fclose(fid);
  end
end


% now lets repeat that, but this time merge all the PDFs
pdf_list = zeros([npmax,maxbin]);
bins_list = pdf_list;

zerobin = 1 + (maxbin-1)/2;            % location of zero bin

disp('summing into 1 PDF for each k...')
for t = times
  tstr=sprintf('%10.4f',t+10000);
  fname=[name,tstr(2:10),'.cpdf'];

  fid=fopen(fname,'r','l');
  if (fid>0) 

    time=fread(fid,1,'float64');
    npmax=fread(fid,1,'float64');         
    for k=1:npmax
      [n_del,delta,bin_size,n_bin,n_call,bins,pdf]=read1pdf(fid);
      % normalize so sum(bin_size * pdf ) = 1
      % this is needed so we can plot PDFs with different binsizes on same axis:
      pdf = pdf / bin_size;
      
      % example: merged pdf has 101 bins  (zero is bin 50)
      % new pdf has 11 bins (zero is bin 5)
      % pdf_list(k,45:55)) = pdf(:)
      x1=zerobin-(n_bin-1)/2;
      x2=zerobin+(n_bin-1)/2;
      pdf_list(k,x1:x2) = pdf_list(k,x1:x2) + pdf(:)'/count;
      bins_list(k,x1:x2) = bins(:);
    end
    fclose(fid);
  end
end

%
% now lets plot the merged pdf:
%
disp('plotting...')
figure(3); clf; 
for k=1:npmax
  pdf = pdf_list(k,:)';
  bin_size = bin_size_list(k);
  bins = bins_list(k,:);
  
  c1=sum(bin_size * pdf.*bins);
  c2=sum(bin_size * pdf.*(bins-c1).^2);
  c4=sum(bin_size * pdf.*(bins-c1).^4);
  flatness(k)=c4/(c2^2);
  
  figure(3);
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
figure(3); hold off;

figure(4); clf;
plot(1:npmax,flatness(1:npmax),'r*') ; hold on
hold off
title('flatness delta filtered')
ax=axis; axis([ax(1),ax(2),1,4]);


%orient tall
figure(3)
print('-dpng',['deltapdf','.png']); 
figure(4)
print('-dpng',['flatness','.png']); 

