%
% read local gradu and gradu**2 ascii files produced by
% extract_gradu.c 
%
% The ascii file contains 8x8x8 subcubes.  For each subcube,
% we have the gradu(3,3) matrix.
%
%

dir='/ccs/scratch/taylorm/dns/decay/';
%dir='/home/mt/';
filename='decay2048-new.0000.7019';
%filename='decay2048-new.0000.7536';

gradu=zeros([3,3,8,8,8]);
gradu2=zeros([3,3,8,8,8]);
fname=[dir,filename,'.gradu'];
fid=fopen(fname);

while 1
  [data,count]=fscanf(fid,'%d',[1,3]);
  if (count==0) break; end;
  i=1+data(1);j=1+data(2);k=1+data(3);
  
  gradu(:,:,i,j,k)=fscanf(fid,'%f',[3,3]);
  
end
fclose(fid);


fname=[dir,filename,'.gradu2'];
fid=fopen(fname);
while 1
  [data,count]=fscanf(fid,'%d',[1,3]);
  if (count==0) break; end;
  i=1+data(1);j=1+data(2);k=1+data(3);
  gradu2(:,:,i,j,k)=fscanf(fid,'%f',[3,3]);
end
fclose(fid);



len=size(gradu2);
epsilon=sum(sum(sum(sum(sum(gradu2)))));
epsilon=epsilon/len(3)/len(4)/len(5);
mu=3.4424e-6;
epsilon=mu*epsilon;
% .7018 time: 2048^3 epsilon=.035
epsilon_l=mu*squeeze(sum(sum(gradu2,1),2));

% normalize:
tscale=(epsilon^(1/3));
gradu=gradu/tscale;
%gradu2=gradu2/tscale/tscale;


mn = min(min(min(min(min(gradu)))));
mx = max(max(max(max(max(gradu)))));
mx=max(abs(mx),abs(mn));
mn=-mx;
mn=.9*mn;
mx=.9*mx;

diag1=reshape(gradu(1,1,:,:,:),[8*8*8,1]);
diag2=reshape(gradu(2,2,:,:,:),[8*8*8,1]);
diag3=reshape(gradu(3,3,:,:,:),[8*8*8,1]);
diag1=[diag1;diag2;diag3];

diag_std=std(diag1);
small=1.4*diag_std;
large=1.90.*small;




%diag_std=1;
figure(2); subplot(2,2,1)
h=hist((diag1/diag_std),25);
hist((diag1/diag_std),25)
title('PDF of <u_{1,1}>')
xlabel('standard deviations u_{1,1}');
hold on;
x=[large,large]/diag_std;
y=[0,.9*max(h)];
plot(x,y,'r',-x,y,'r');
x=[small,small]/diag_std;
plot(x,y,'g',-x,y,'g');
hold off;


diag1=reshape(gradu(1,2,:,:,:),[8*8*8,1]);
diag2=reshape(gradu(1,3,:,:,:),[8*8*8,1]);
diag3=reshape(gradu(2,1,:,:,:),[8*8*8,1]);
diag4=reshape(gradu(2,3,:,:,:),[8*8*8,1]);
diag5=reshape(gradu(3,1,:,:,:),[8*8*8,1]);
diag6=reshape(gradu(3,2,:,:,:),[8*8*8,1]);
offdiag1=[diag1;diag2;diag3;diag4;diag5;diag6];
figure(2); subplot(2,2,3)
h=hist((diag1/diag_std),25);
hist((diag1/diag_std),25);
title('PDF of <u_{1,2}>')
xlabel('standard deviations u_{1,1}');
hold on;
x=[large,large]/diag_std;
y=[0,.9*max(h)];
plot(x,y,'r',-x,y,'r');
x=[small,small]/diag_std;
plot(x,y,'g',-x,y,'g');
hold off;

figure(2); subplot(2,2,4)
diag1=reshape(epsilon_l,[8*8*8,1]);
h=hist(diag1,25);
hist(diag1,25);
title('PDF local \epsilon')
hold on;
x=[epsilon,epsilon];
y=[0,.9*max(h)];
plot(x,y,'r')
hold off;


figure(2); subplot(2,2,2)
image(squeeze(gradu(:,:,1,1,1)),'CDataMapping','scaled' );
title('\nabla u \epsilon^{-1/3} L^{2/3}')
caxis([mn,mx])
colorbar
print -dpsc -r72  gradu-stats.ps
print -djpeg -r72  gradu-stats.jpg






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
%      if (mxa/mx_remaining > 2.40) 
      if (mxa>large & mx_remaining < small) 
        count=count+1;
        if (count==1) clf; end;
        subplot(m1,m2,count);
        image(squeeze(gradu(:,:,i,j,k)),'CDataMapping','scaled' );
        caxis([mn,mx])
        set(gca,'YTickLabel','')
        set(gca,'XTickLabel','')
        axis image
        xlabel(sprintf('(%i,%i,%i)',i-1,j-1,k-1))

        if (count==m1*m2)
          count=0;
          disp('warning: more than 25 subcubes selected. only the last plot will be saved');
        end


      end 
    end
  end
end

subplot(m1,m2,1)
text(.5,0,sprintf('threshold(std)=%.2f, %.2f',small/diag_std,large/diag_std));
figure(1);
print('-dpsc',sprintf('gradu-set.ps'));
print('-djpeg','-r72',sprintf('gradu-set.jpg'));




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
      

      subplot(m1,m2,count);
      if (count==1) clf; end;
      image(squeeze(gradu(:,:,i,j,k)),'CDataMapping','scaled' );
      caxis([mn,mx])
      set(gca,'YTickLabel','')
      set(gca,'XTickLabel','')
      axis image

      if (count==m1*m2)
        count=0;
        figcount=figcount+1
        print('-dpsc',sprintf('gradu-%i.ps',figcount)); 
        print('-djpeg','-r72',sprintf('gradu-%i.jpg',figcount)); 
      end

    end
  end
end
