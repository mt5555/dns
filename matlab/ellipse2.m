%
%########################################################################
%#  plotting ellipses
%########################################################################

%ts=input('time=? ');

%name='../src/temp';
%name='../src/temp0000.0000.ellipse';
%name='../src/vxpair/vx4096b0009.0000.ellipse';
%name='/home/taylorm/vxpair/vx2048a0050.0000.ellipse';
%name='/data/vxpair/vx2048c0000.0000.ellipse';

name='/ccs/taylorm/dns/src/vxpair/vx6144c';
times=[0:.1:.5];


for t=times
  tstr=sprintf('%10.4f',t+10000);
  fname=[name,tstr(2:10),'.ellipse2'];
  disp(fname)
  fid=endianopen(fname,'r');
  if (fid>=0) 

     [nell,count]=fread(fid,1,'float64');
     if (count~=1) break; end;
     np=fread(fid,1,'float64');
     time=fread(fid,1,'float64');

     npv=(1:np)';
     cosc = cos(2*pi*(npv-1)/(np));
     sinc = sin(2*pi*(npv-1)/(np));
     cos2c = cos(2*pi*2*(npv-1)/(np));
     sin2c = sin(2*pi*2*(npv-1)/(np));


     figure(1); clf; 

     disp('nell      Rmin          Rmax        m=1/m0       m=2/m0');

     for i=1:nell
        center=fread(fid,2,'float64');
        plot(center(1),center(2),'.'); hold on; 
        [rad,count]=fread(fid,np,'float64');
        x=center(1) + rad.*cosc;
        y=center(2) + rad.*sinc;
        x(np+1)=x(1);
        y(np+1)=y(1);
        plot(x,y)
        %plot(x,y,'r.')

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

      
     end     
     axis([1 3 0 1.5]);
     axis equal
     title(sprintf('time=%f',time))
     hold off;
     %'pause' ;     pause
     fclose(fid);
  end
end

return


