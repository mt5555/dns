%
% read local gradu and gradu**2 ascii files produced by
% extract_gradu.c 
%
% The ascii file contains 8x8x8 subcubes.  For each subcube,
% we have the gradu(3,3) matrix.
%
%

dir='/ccs/scratch/taylorm/dns/decay/';
filename='decay2048-new.0000.7019';
%filename='decay2048-new.0000.7536';

gradu=zeros([3,3,8,8,8]);
fname=[dir,filename,'.gradu'];
fid=fopen(fname);

while 1
  [data,count]=fscanf(fid,'%d',[1,3]);
  if (count==0) break; end;
  i=1+data(1);j=1+data(2);k=1+data(3);
  
  gradu(:,:,i,j,k)=fscanf(fid,'%f',[3,3]);
  
end
fclose(fid);

mn = min(min(min(min(min(gradu)))));
mx = max(max(max(max(max(gradu)))));
mx=max(abs(mx),abs(mn));
mn=-mx;

diag1=reshape(gradu(1,1,:,:,:),[8*8*8,1]);
diag2=reshape(gradu(2,2,:,:,:),[8*8*8,1]);
diag3=reshape(gradu(3,3,:,:,:),[8*8*8,1]);
diag1=[diag1;diag2;diag3];
diag_std=std(diag1);
small=1.5*diag_std;
large=2*diag_std;







figure(2); subplot(1,1,1)
image(squeeze(gradu(:,:,1,1,1)),'CDataMapping','scaled' );
caxis([mn,mx])
colorbar
print -dpsc -r72  gradu-colormap.ps
print -djpeg -r72  gradu-colormap.jpg



m1=5;
m2=5;
figure(1); clf;
figcount=0;
count=0;
for i=1:8
  for j=1:8
    for k=1:8
      a=squeeze(gradu(:,:,i,j,k));
      offdiag=a-diag(diag(a));

      [a1,mxi]=max(offdiag); [mxa,mxj]=max(a1); mxi=mxi(mxj);
      [a1,mni]=min(offdiag); [mna,mnj]=min(a1); mni=mni(mnj);

      if (mxa < -mna ) 
        mxa=mna;
        mxi=mni;
        mxj=mnj;
      end
      a(mxi,mxj)=0;
      mx_remaining=max(max(abs(a)));
      mxa=abs(mxa);
      if (mxa/mx_remaining > 2.40) 
        count=count+1;
        if (count>m1*m2)
          count=1;
          figcount=figcount+1
          print('-dpsc',sprintf('gradu-%i.ps',figcount)); 
          print('-djpeg','-r72',sprintf('gradu-%i.jpg',figcount)); 
        end
        
        subplot(m1,m2,count);
        image(squeeze(gradu(:,:,i,j,k)),'CDataMapping','scaled' );
        caxis([mn,mx])
        set(gca,'YTickLabel','')
        set(gca,'XTickLabel','')
        axis image
        xlabel(sprintf('(%i,%i,%i)',i-1,j-1,k-1))
      end 
    end
  end
end

return



m1=8;
m2=8;
figure(3);
figcount=0;
count=0;
for i=1:8
  for j=1:8
    for k=1:8
      count=count+1;
      
      if (count>m1*m2)
        count=1;
        figcount=figcount+1
        print('-dpsc',sprintf('gradu-%i.ps',figcount)); 
        print('-djpeg','-r72',sprintf('gradu-%i.jpg',figcount)); 
      end


      subplot(m1,m2,count);
      image(squeeze(gradu(:,:,i,j,k)),'CDataMapping','scaled' );
      caxis([mn,mx])
      set(gca,'YTickLabel','')
      set(gca,'XTickLabel','')
      axis image

    end
  end
end

