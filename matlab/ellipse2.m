%
%########################################################################
%#  plotting ellipses
%########################################################################
clear;

%ts=input('time=? ');

%name='../src/temp';
%name='../src/temp0000.0000.ellipse';
%name='../src/vxpair/vx4096b0009.0000.ellipse';
%name='/home/taylorm/vxpair/vx2048a0050.0000.ellipse';
%name='/data/vxpair/vx2048c0000.0000.ellipse';

name='/ccs/taylorm/dns/src/vxpair/vx6144e';
times=[17.28:.02:30.0];

%name='/ccs/taylorm/dns/src/vxpair/vx6144c';
%times=[10.00:.1:30.0];

%name='/ccs/taylorm/dns/src/vxpair/vx6144d';
%times=[0.00:.02:50.0];

ccol=[ 'b','g','r','c','m','y', 'b','g','r','c','m','y' ];  

emode=[];
timev=[];
k=0;
for t=times
  tstr=sprintf('%10.4f',t+10000);
  fname=[name,tstr(2:10),'.ellipse2'];
  fid=endianopen(fname,'r');


  
  
  if (fid>=0) 
     disp(fname)

     [nell,count]=fread(fid,1,'float64');
     if (count~=1) break; end;
     np=fread(fid,1,'float64');
     time=fread(fid,1,'float64');
     wcenter=fread(fid,2,'float64');
     
     npv=(1:np)';
     cosc = cos(2*pi*(npv-1)/(np));
     sinc = sin(2*pi*(npv-1)/(np));
     cos2c = cos(2*pi*2*(npv-1)/(np));
     sin2c = sin(2*pi*2*(npv-1)/(np));


     figure(1);  subplot(1,1,1) ; hold on;
     clf; 
     plot(wcenter(1),wcenter(2),'r.'); hold on;
     
     disp('nell      Rmin          Rmax        m=1/m0       m=2/m0');

     k=k+1;
     timev=[timev,time];
     for i=1:nell
        wval(i)=fread(fid,1,'float64');
        center=fread(fid,2,'float64');
        plot(center(1),center(2),'.'); hold on; 
        [rad,count]=fread(fid,np,'float64');

        npnew=np;

        %[i,np] 
        if (np==32000) 
           rad=rad(1:2:np);

           npnew=length(rad);
           npv=(1:npnew)';
           cosc = cos(2*pi*(npv-1)/(npnew));
           sinc = sin(2*pi*(npv-1)/(npnew));
           cos2c = cos(2*pi*2*(npv-1)/(npnew));
           sin2c = sin(2*pi*2*(npv-1)/(npnew));
        end            

        x=center(1) + rad.*cosc;
        y=center(2) + rad.*sinc;
        
        x(npnew+1)=x(1);
        y(npnew+1)=y(1);
        %plot(x,y)
        plot(x,y,[ccol(i),'.'])

        % compute FFT
        sq2=sqrt(2d0);

        dft(1)=sum(rad);
        dft(2)=sum(rad.*cosc*sq2);
        dft(3)=sum(rad.*sinc*sq2);
        dft(4)=sum(rad.*cos2c*sq2);
        dft(5)=sum(rad.*sin2c*sq2);

        disp(sprintf('%i      %f     %f      %f     %f ',i,min(rad),max(rad),...
             sqrt(dft(2)^2+dft(3)^2)/dft(1),...
             sqrt(dft(4)^2+dft(5)^2)/dft(1)))

        if (length(emode)==0) 
           emode=zeros([nell,2,1]);
        end
        emode(i,1,k)=sqrt(dft(4)^2+dft(5)^2);
        emode(i,2,k)=sqrt(dft(2)^2+dft(3)^2);
     end     
     fclose(fid);

     % try and get some vorticity contours:
     fname=[name,tstr(2:10),'.vor'];
     [x,y,z,vor,time2]=getfield(fname);
     if (time2>0) 
       disp('plotting vorticity contours');
       vor = squeeze(vor(:,:,1));
       subsample=2;
       if (subsample>1) 
         nx=length(x);
         ny=length(y);
         vor=vor(1:subsample:nx,1:subsample:ny);
         x=x(1:subsample:nx);
         y=y(1:subsample:ny);
       end
       v=  [3.3938    2.2626    1.1313    0.5656];
       contour(x,y,vor',v)
     end
     axis equal
     axis([1 3 0 1.5]);
     title(sprintf('time=%f',time))
     hold off;

     figure(2); clf;
     for i=1:nell
       subplot(2,1,1)
       semilogy(timev,squeeze(emode(i,1,:)),[ccol(i),'-']); hold on;
       subplot(2,1,2)
       semilogy(timev,squeeze(emode(i,2,:)),[ccol(i),'-']); hold on;
     end
     hold off;
     'pause' ;     pause
     
  end
end

      

return


