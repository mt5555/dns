function ...
    [nx,ndelta,ndir,r_val,ke,epsilon,mu,...
    D_ll,D_lll,D1_tt,D2_tt,D1_ltt,D2_ltt,...
    SP_lll,SN_lll,SP1_ltt,SP2_ltt,SN1_ltt,SN2_ltt,H_ltt,H_tt,D_lltt,Dl,Dt,...
    h_epsilon,Q_epsilon] ...
     = readisostr(fname)


fid=endianopen(fname,'r')

fname 

ndelta=fread(fid,1,'float64')
ndir  =fread(fid,1,'float64')
nlon  =fread(fid,1,'float64')
ntran =fread(fid,1,'float64')
nscalars =fread(fid,1,'float64')
nnew2 =fread(fid,1,'float64')

SP_lll=0;
SN_lll=0;
SP1_ltt=0;
SN1_ltt=0;
SP2_ltt=0;
SN2_ltt=0;


r_val=fread(fid,[ndelta,ndir],'float64');
if (nlon==2) 
   D_ll=fread(fid,[ndelta,ndir],'float64');
   D_lll=fread(fid,[ndelta,ndir],'float64');
end
if (nlon==4) 
   D_ll=fread(fid,[ndelta,ndir],'float64');
   D_lll=fread(fid,[ndelta,ndir],'float64');
   SP_lll=fread(fid,[ndelta,ndir],'float64');
   SN_lll=fread(fid,[ndelta,ndir],'float64');
end
if (nlon==11) 
   Dl=zeros([ndelta,ndir,9]);
   for p=2:10 
     temp=fread(fid,[ndelta,ndir],'float64');
     Dl(:,:,p-1)=temp;
     if (p==2) D_ll=temp; end;
     if (p==3) D_lll=temp; end;
   end
   SP_lll=fread(fid,[ndelta,ndir],'float64');
   SN_lll=fread(fid,[ndelta,ndir],'float64');
end
if (nlon==12) 
   Dl=zeros([ndelta,ndir,9]);
   for p=2:10 
     temp=fread(fid,[ndelta,ndir],'float64');
     Dl(:,:,p-1)=temp;
     if (p==2) D_ll=temp; end;
     if (p==3) D_lll=temp; end;
   end
   SP_lll=fread(fid,[ndelta,ndir],'float64');
   SN_lll=fread(fid,[ndelta,ndir],'float64');
   H_ltt=fread(fid,[ndelta,ndir],'float64');
end
if (nlon==13) 
   Dl=zeros([ndelta,ndir,9]);
   for p=2:10 
     temp=fread(fid,[ndelta,ndir],'float64');
     Dl(:,:,p-1)=temp;
     if (p==2) D_ll=temp; end;
     if (p==3) D_lll=temp; end;
   end
   SP_lll=fread(fid,[ndelta,ndir],'float64');
   SN_lll=fread(fid,[ndelta,ndir],'float64');
   H_ltt=fread(fid,[ndelta,ndir],'float64');
   H_tt=fread(fid,[ndelta,ndir],'float64');
end
if (nlon==14) 
   Dl=zeros([ndelta,ndir,9]);
   for p=2:10 
     temp=fread(fid,[ndelta,ndir],'float64');
     Dl(:,:,p-1)=temp;
     if (p==2) D_ll=temp; end;
     if (p==3) D_lll=temp; end;
   end
   SP_lll=fread(fid,[ndelta,ndir],'float64');
   SN_lll=fread(fid,[ndelta,ndir],'float64');
   H_ltt=fread(fid,[ndelta,ndir],'float64');
   H_tt=fread(fid,[ndelta,ndir],'float64');
   D_lltt=fread(fid,[ndelta,ndir],'float64');
end

if (ntran==2) 
   D1_tt=fread(fid,[ndelta,ndir],'float64');
   D2_tt=fread(fid,[ndelta,ndir],'float64');
   D1_ltt=fread(fid,[ndelta,ndir],'float64');
   D2_ltt=fread(fid,[ndelta,ndir],'float64');
end
if (ntran==4) 
   D1_tt=fread(fid,[ndelta,ndir],'float64');
   D2_tt=fread(fid,[ndelta,ndir],'float64');
   D1_ltt=fread(fid,[ndelta,ndir],'float64');
   D2_ltt=fread(fid,[ndelta,ndir],'float64');
   SP1_ltt=fread(fid,[ndelta,ndir],'float64');
   SP2_ltt=fread(fid,[ndelta,ndir],'float64');
   SN1_ltt=fread(fid,[ndelta,ndir],'float64');
   SN2_ltt=fread(fid,[ndelta,ndir],'float64');
end
if (ntran==12) 
  Dt1=zeros([ndelta,ndir,9]);
  Dt2=zeros([ndelta,ndir,9]);
  for p=2:10
   temp1=fread(fid,[ndelta,ndir],'float64');    
   temp2=fread(fid,[ndelta,ndir],'float64');    
   Dt1(:,:,p-1)=temp1;
   Dt2(:,:,p-1)=temp2;
   if (p==2)
     D1_tt=temp1;
     D2_tt=temp2;
   end
 end 
   D1_ltt=fread(fid,[ndelta,ndir],'float64');
   D2_ltt=fread(fid,[ndelta,ndir],'float64');
   SP1_ltt=fread(fid,[ndelta,ndir],'float64');
   SP2_ltt=fread(fid,[ndelta,ndir],'float64');
   SN1_ltt=fread(fid,[ndelta,ndir],'float64');
   SN2_ltt=fread(fid,[ndelta,ndir],'float64');
end


if (nscalars>=7) 
  time=fread(fid,1,'float64');
  nx=fread(fid,1,'float64');
  ny=fread(fid,1,'float64');
  nz=fread(fid,1,'float64');
  mu=fread(fid,1,'float64');   
  ke=fread(fid,1,'float64');
  epsilon=fread(fid,1,'float64'); 

  h_epsilon=1;
  if (nscalars>=8) 
    [h_epsilon,count]=fread(fid,1,'float64');  
  end
  
  Q_eps=1;
  if (nscalars>=9)
     [Q_eps,count]=fread(fid,1,'float64');
  end

  eta = (mu^3 / epsilon)^.25;
  delx_over_eta=(1/nx)/eta;
end

tmp = fread(fid,1,'float64');
tmp=size(tmp);
if (tmp(1)~=0) 
  disp('Error reading input file...')
end
fclose(fid);



% D_lll=SP_lll-SN_lll;
% but to make SP and SN have the same sign as D, multipl
% by -1 so:  D=SN-SP
SP_lll=-SP_lll;
SN_lll=-SN_lll;
SP1_ltt=-SP1_ltt;
SN1_ltt=-SN1_ltt;
SP2_ltt=-SP2_ltt;
SN2_ltt=-SN2_ltt;

Dt=.5*(Dt1+Dt2);
