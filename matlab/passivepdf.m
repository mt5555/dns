%
%########################################################################
%  read and plot passive scalar pdf's
%########################################################################
%
%
clear all;


fid=fopen('/home/taylorm/data/dns/tmix/tmix256D/tmix256D0001.0000.spdf','r','l');
%fid=fopen('../src/temp0000.0000.spdf','r','l');

time=fread(fid,1,'float64');
npmax=fread(fid,1,'float64');         
disp(sprintf('npassive = %i',npmax))

figure(1); clf; subplot(npmax/2,2,1)

np=1;
for p=1:npmax
   [n_del,delta,bin_size,n_bin,n_call,bins,pdf]=read1pdf(fid);
   %if (p==np) break; end;

   c1=sum(pdf.*bins);
   c2=sum(pdf.*(bins-c1).^2);
   c4=sum(pdf.*(bins-c1).^4);

%   c2=sum(pdf.*(bins).^2);
%   c2=c2-c1.^2
   
   
   mx=max(bins - bins.*(pdf==0));
   mn=min(bins - bins.*(pdf==0));        % min over non zero values
   
   subplot(npmax/2,2,p) 
   plot(bins,pdf)
   ax=axis;
   axis([-.1, 1.1, 0, .15]);
   xlabel(sprintf('<c^2>=%.4f  K=%.4f ',c2,c4/c2/c2));
   set(gca,'YTickLabel','')    
   if (p==1) 
      title(sprintf('t=%.4f',time)); 
   end

   

end
times=sprintf('%.4f',time+10000);
times=times(2:length(times));
orient tall
print('-dpsc',['ppdf',times,'.ps']); 
%print('-djpeg',['ppdf',times,'.jpg']); 
