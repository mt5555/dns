%
%########################################################################
%#  plot passive scalar
%########################################################################
%

name = '/scratch2/taylorm/tmix256C/tmix256C'
time=1.10;
times=sprintf('%.5f',time+10000);
times=times(2:length(times)-1);

schmidt_list=[.01 .05 .1 .5 1];
type_list=0:1;
npassive=length(schmidt_list)*length(type_list);

k=0;
for sch=schmidt_list
for type=type_list

  ext=sprintf('%3i',type+100);
  ext=ext(2:length(ext));
  ext2=sprintf('%8.3f',sch+1000);
  ext2=ext2(2:length(ext2));

  ext=['.t',ext,'.s',ext2];
  
  s=findstr(name,'/');
  s=s(length(s));
  shortname=name(s+1:length(name));

  k=k+1;

  % try to read the pscalars-turb file 
  fname=[name,times,'.pscalars-turb'];
  disp(fname)
  lambda_c=0;
  fid=endianopen(fname,'r');
  if (fid>=0) 
    [ns_e,count] = fread(fid,1,'float64');
    [npassive,count] = fread(fid,1,'float64');
    time_2 = fread(fid,1,'float64');
    mu = fread(fid,1,'float64');
    pints_e=zeros([2+ns_e,npassive,1]); 
    for np=1:npassive
      data1 = fread(fid,[ns_e,1],'float64');
      data1=[time;mu;data1];
      pints_e(:,np)= data1;
    end
    fclose(fid);
    np=k; 
    c1=squeeze(pints_e(26,np))';
    c2=squeeze(pints_e(4,np))';        % index=2 
    c2=c2-c1.^2; 
    cx2=zeros([3,length(c1)]);
    cx2(1,:)=pints_e(5,np);
    cx2(2,:)=pints_e(6,np);
    cx2(3,:)=pints_e(7,np);
    lambda_c=sqrt(c2./mean(cx2))

% $$$     % try to read the scalars-turb file 
% $$$     fname=[name,times,'.scalars-turb'];
% $$$     disp(fname)
% $$$     fid=endianopen(fname,'r');
% $$$     [ns_e,count] = fread(fid,1,'float64');
% $$$     time_2 = fread(fid,1,'float64');
% $$$     data1 = fread(fid,[ns_e,1],'float64');
% $$$     data1=[time;data1];
% $$$     ints_e= data1;
% $$$     fclose(fid);
% $$$     for i=1:3
% $$$       ux2(i)=ints_e(i+1);    % < u1,1 ^ 2 >
% $$$     end
% $$$     epsilon=15*mu.*mean(ux2);
% $$$     eta = (mu.^3./epsilon).^(.25);
% $$$     eta_c= eta/sqrt(sch)
  end

  
  
  
  
   fname=[name,times,ext]
   [x,y,z,s,time]=getfield(fname);
   if (time<0) 
     disp('error opening file ');
     return;
   end
   smax=max(max(s));
   [mx,slice1]=max(smax);
   smax=min(min(s));
   [mn,slice2]=min(smax);

   figure(1);
   subplot(npassive/2,2,k)
   splot=squeeze(s(:,:,slice1));
   pcolor(x,y,splot')
   
   stitle=sprintf('%s    time=%.2f  max=%f',shortname,time,mx)
   if (k==1) title(stitle); end;
   axis equal
   axis([0,max(x),0,max(y)]);
   shading interp
   caxis([0 1]) 

   

   if (lambda_c>0) 
     % black box of size 5 lambda_c
     fac=2*pi; facl='2 \pi \lambda_c';
     len=fac*lambda_c;
     xt=[1,1-len,1-len,1,1];
     yt=[0,0,len,len,0];
     text(1.05,yt(3)/2,facl);
     %patch(xt,yt,'w')
     line(xt,yt,'LineWidth',1.5,'Color','k');

     figure(2); 
     subplot(npassive/2,2,k)
     splot=squeeze(s(:,:,slice1));
     pcolor(x,y,splot')
     
     stitle=sprintf('%s    time=%.2f  max=%f',shortname,time,mx)
     if (k==1) title(stitle); end;
     axis equal
     axis([0,fac*lambda_c,0,fac*lambda_c]);
     shading interp
     caxis([0 1]) 
   end
   
   
   
end
end
figure(1)
orient tall
print('-djpeg','-r125',['p',times,'.jpg']); 
figure(2)
orient tall
print('-djpeg','-r125',['pz',times,'.jpg']); 

