%
%########################################################################
%#  plot passive scalar
%########################################################################
%
figure(1);clf;
clear
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

  k=k+1;
  ext=sprintf('%3i',type+100);
  ext=ext(2:length(ext));
  ext2=sprintf('%8.3f',sch+1000);
  ext2=ext2(2:length(ext2));

  ext=['.t',ext,'.s',ext2];
  
  s=findstr(name,'/');
  s=s(length(s));
  shortname=name(s+1:length(name));

  
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
    np=11-k; 
    if (sch ~= pints_e(3,np,1))
      [sch,pints_e(3,np,1)] 
      disp('wrong schmidt number in scalars file')
      return 
    end  
    c1=squeeze(pints_e(26,np))';
    c2=squeeze(pints_e(4,np))';        % index=2 
    c2=c2-c1.^2; 
    cx2=zeros([3,length(c1)]);
    cx2(1,:)=pints_e(5,np);
    cx2(2,:)=pints_e(6,np);
    cx2(3,:)=pints_e(7,np);
    lambda_c=sqrt(c2./mean(cx2))

    % try to read the scalars-turb file 
    fname=[name,times,'.scalars-turb'];
    disp(fname)
    fid=endianopen(fname,'r');
    [ns_e,count] = fread(fid,1,'float64');
    time_2 = fread(fid,1,'float64');
    data1 = fread(fid,[ns_e,1],'float64');
    data1=[time;data1];
    ints_e= data1;
    fclose(fid);
    for i=1:3
      ux2(i)=ints_e(i+1);    % < u1,1 ^ 2 >
    end
    epsilon=15*mu.*mean(ux2);
    eta = (mu.^3./epsilon).^(.25);
    eta_c= eta/sqrt(sch)
  end

  
  
  
   fname=[name,'-gradxy2',times,ext]
   [x,y,z,s,time]=getfield(fname);
   smax=max(max(s));
   [mx,slice1]=max(smax);
   smax=min(min(s));
   [mn,slice2]=min(smax);

   splot=squeeze(s(:,:,slice1));

   fac=2*pi; facl='2 \pi \eta_c';
   len=min(1,fac*eta_c);

   figure(2)
   subplot(npassive/2,2,k)
   pcolor(x,y,log(splot'));  caxis([log(.01) log(13)]) ;
   axis equal
   axis([0,max(x),0,max(y)]);
   shading interp


   % patch of size 10eta_c
   xt=[1,1-len,1-len,1,1];
   yt=[0,0,len,len,0];
   patch(xt,yt,'w');
   text(1.05,yt(3)/2,facl);

   
   figure(4)
   subplot(npassive/2,2,k)
   pcolor(x,y,log(splot'));  caxis([log(.01) log(13)]) ;
   axis equal
   axis([0,len,0,len]);
   shading interp
   

   len=min(1,fac*eta_c);
if (0)   
   figure(1)
   subplot(npassive/2,2,k)
   pcolor(x,y,splot');   caxis([0 13]) ;
   axis equal
   axis([0,max(x),0,max(y)]);
   shading interp
   xt=[1,1-len,1-len,1,1];
   yt=[0,0,len,len,0];
   patch(xt,yt,'w');
   text(1.05,yt(3)/2,facl);

   stitle=sprintf('time=%.2f  min=%f  max=%f',time,mn,mx)
   if (k==1) title(stitle); end;


   figure(3)
   subplot(npassive/2,2,k)
   pcolor(x,y,splot');   caxis([0 13]) ;
   axis equal
   axis([0,len,0,len]);
   shading interp
end      
   

end
end
%figure(1)
%orient tall
%print('-djpeg','-r125',['pdiss',times,'.jpg']); 
%figure(3)
%orient tall
%print('-djpeg','-r125',['pdiss_z',times,'.jpg']); 

figure(2)
orient tall
print('-djpeg','-r125',['plogdiss',times,'.jpg']); 

figure(4)
orient tall
print('-djpeg','-r125',['plogdiss_z',times,'.jpg']); 


