#include "macros.h"
module ellipse
use params
implicit none
#if 0

module which will compute ellipse centerered around the peak vorticity


#endif



integer,private   :: init = 0

integer,parameter :: nelld = 7  !  number of ellipses
integer,parameter :: npd   = 65  !  number of points along each ellipse


real*8 :: wval(nelld)         ! vorticity contour values (% of max vorticity)
real*8 :: cosc(npd),sinc(npd)     ! angles to perform contouring
real*8 :: cos2c(npd),sin2c(npd)   
real*8 :: Rad(npd,nelld)       ! radius
real*8 :: center(2)            ! location of center
real*8 :: Rad2(npd,nelld)       ! shifted radius
real*8 :: ccord(2,nelld)        ! shifted centers
real*8 :: mxw                 ! max vorticity on grid
real*8 :: mxw_init=-1         ! max vorticity at time=0
real*8 :: dft(0:4,nelld)              ! modes of Rad

real*8 :: contour_eps = 1e-6     ! find contours to within this accuracy
real*8 :: center_eps  = 1e-5    ! find center to within this accuracy




contains

subroutine ellipse_init
implicit none
integer :: nell,np
init=1


wval(1)=7/8.
wval(2)=6/8.
wval(3)=5/8.
wval(4)=4/8.
wval(5)=3/8.
wval(6)=2/8.
wval(7)=1/8.
if (nelld /= 7) then
   call abort("ellipse init error")
endif

do np=1,npd
   cosc(np) = cos(2*pi*(np-1)/(npd))
   sinc(np) = sin(2*pi*(np-1)/(npd))
   cos2c(np) = cos(2*2*pi*(np-1)/npd)
   sin2c(np) = sin(2*2*pi*(np-1)/npd)
enddo
Rad=0


end subroutine






subroutine ellipse_output(time)
use params
implicit none
real*8 :: xell(npd),yell(npd)
real*8 :: tmp,time
integer :: nell,np,ierr
CPOINTER fid
character(len=280) :: fname
character(len=80) :: message


if (io_pe==my_pe) then
   print *,'center: ',center
   print *,'nell      Rmin          Rmax        m=1/m0        m=2/m0'
   do nell=1,nelld
      write(*,'(i1,2f14.8,f14.8,f14.8)') nell,minval(Rad2(:,nell)),maxval(Rad2(:,nell)),&
          sqrt(dft(1,nell)**2+dft(2,nell)**2)/dft(0,nell),&
          sqrt(dft(3,nell)**2+dft(4,nell)**2)/dft(0,nell)
   enddo

   write(message,'(f10.4)') 10000.0000 + time
   fname = rundir(1:len_trim(rundir)) // runname(1:len_trim(runname)) // message(2:10) // ".ellipse2"
   call copen(fname,"w",fid,ierr)
   if (ierr/=0) then
      write(message,'(a,i5)') "spec_write(): Error opening file errno=",ierr
      call abort(message)
   endif

   
   tmp=nelld
   call cwrite8(fid,tmp,1)
   tmp=npd
   call cwrite8(fid,tmp,1)
   call cwrite8(fid,time,1)

   do nell=1,nelld
!      do np=1,npd
!         xell(np)=center(1)+Rad(np,nell)*cosc(np)
!         yell(np)=center(2)+Rad(np,nell)*sinc(np)
!      enddo
!      call cwrite8(fid,xell,npd)
!      call cwrite8(fid,yell,npd)
      call cwrite8(fid,ccord(1,nell),2)
      call cwrite8(fid,Rad2(1,nell),npd)
   enddo
   call cclose(fid,ierr)
endif

end subroutine







subroutine comp_ellipse(w,setmax)
use params
use mpi
implicit none
real*8 :: w(nx,ny)
integer :: setmax

!local
integer :: nell,np
real*8 :: tmx1,tmx2

call wallclock(tmx1)
if (init==0) call ellipse_init()

! estimate center
call findcenter(w,center)
if (setmax==1) then
   mxw_init=mxw
   return
endif

if (mxw_init<=0) then
   call abort("ellipse(): mxw_init was not set")
endif



! compute coordinates of ellipse:  mxcord + Rad exp(i 2pi theta)
! uses binary search and 4th order interpolation
do nell=1,nelld
   call findellipse(w,nell,center,Rad(1,nell))
enddo

! iterate to find the "true" center of the ellipse
call findbestcenter(w,center)


call wallclock(tmx2)
tims(17)=tims(17)+(tmx2-tmx1)

end subroutine









subroutine findbestcenter(w,mxcord)
use params
use mpi
implicit none

real*8 :: w(nx,ny),mxcord(2)
real*8 :: xell(npd),yell(npd)
real*8 :: sq2,relax
integer np,nell,count

sq2=sqrt(2d0)
relax=1.50

do nell=1,nelld
   if (wval(nell)*mxw_init<mxw) then

   count=0
50 continue

   Rad2(:,nell)=Rad(:,nell)
   ccord(:,nell)=mxcord
   do np=1,npd
      xell(np)=mxcord(1)+Rad(np,nell)*cosc(np)
      yell(np)=mxcord(2)+Rad(np,nell)*sinc(np)
   enddo


   ! find best mxcord(:) to minimize non-ellipticial modes in Rad
   do
      dft(:,nell)=0
      do np=1,npd
         dft(0,nell)=dft(0,nell) + Rad2(np,nell)
         dft(1,nell)=dft(1,nell) + Rad2(np,nell)*cosc(np)*sq2
         dft(2,nell)=dft(2,nell) + Rad2(np,nell)*sinc(np)*sq2
         dft(3,nell)=dft(3,nell) + Rad2(np,nell)*cos2c(np)*sq2
         dft(4,nell)=dft(4,nell) + Rad2(np,nell)*sin2c(np)*sq2
      enddo
      dft(:,nell)=dft(:,nell)/npd

      if ( (dft(1,nell)**2 + dft(2,nell)**2) < dft(0,nell)*center_eps**2 ) exit
      if (count>40) then
         exit
      else if (count==20) then
         call print_message("findbestcenter(): restarting iteration")
         relax=.9
         count=count+1
         goto 50
      endif

      ! move the center mxcord(:)
      !print *,nell
      !print *,ccord(1),dft(1,nell)
      !print *,ccord(2),dft(2,nell)
      ccord(1,nell)=ccord(1,nell)+relax*dft(1,nell)     ! 3.0 failes, 2.5 works
      ccord(2,nell)=ccord(2,nell)+relax*dft(2,nell)
      count=count+1


      call findellipse(w,nell,ccord(1,nell),Rad2(1,nell))
   enddo
   endif
enddo
end subroutine









subroutine findcenter(w,mxcord)
use params
use mpi
implicit none
real*8 :: w(nx,ny),mxcord(2)

!local
real*8 :: tmp1,tmp2,mxcord2(2)
integer :: i,j,ierr

!
! find max vorticity
!
mxw = -9d20
do j=inty1,inty2
do i=intx1,intx2
   if (w(i,j)>mxw) then
      mxcord(1)=xcord(i)
      mxcord(2)=ycord(j)
      mxw=w(i,j)
   endif
enddo
enddo

#ifdef USE_MPI
   tmp1=mxw
   call MPI_allreduce(tmp1,tmp2,1,MPI_REAL8,MPI_MAX,comm_3d,ierr)
   if (tmp2==mxw) then
      ! we have the maximum.  leave mxcord unchanged
   else
      mxw=tmp2
      mxcord=-9d20
   endif
   mxcord2=mxcord
   call MPI_allreduce(mxcord2,mxcord,2,MPI_REAL8,MPI_MAX,comm_3d,ierr)
#endif
end subroutine






#if 1
subroutine findellipse(w,nell,mxcord,r)
use params
use mpi
implicit none
real*8 :: w(nx,ny),mxcord(2),r(npd)
integer :: nell

!local
real*8 :: mxcord2(2)
real*8 :: wcontour,winterp
real*8 :: Rdelta,tmp(2),tmpout(2)
integer :: np,count,ierr

!
! compute location of ellipse:  mxcord(1) + r(np) * cosc(np) 
!                               mxcord(2) + r(np) * sinc(np) 
!

   wcontour = wval(nell)*mxw_init
   if (wcontour<mxw) then
   !print *,'looking for ',wcontour,mxw
   do np=1,npd
      ! find Rad so that:  w(Rad cosc(np), Rad sinc(np)) = wcontour

      Rdelta=.01
      r(np)=Rdelta

      count=0
      do 
         call interp4w(w,mxcord(1),mxcord(2),r(np),cosc(np),sinc(np),winterp)
         if (winterp<-9d10) then
            r(np)=-9d20  ! not on this CPU
         else if (winterp>wcontour) then
            r(np)=r(np)+Rdelta ! undershoot
         else if (winterp<=wcontour) then ! overshoot
            Rdelta=Rdelta/2
            r(np)=r(np)-Rdelta
         endif
#ifdef USE_MPI
         ! r(np):  one cpu has valid (positive) value, the rest
         ! have -9d20.  so take MAX to update all r(np) on all cpus.
         ! For Rdelta, we want to give all processors the value 
         ! from the cpu which contained the interpolation point. 
         ! Since Rdelta is fixed, or decreasing, taking the MIN will
         ! achive this: 
         tmp(1)=r(np)
         tmp(2)=-Rdelta
         call MPI_allreduce(tmp,tmpout,2,MPI_REAL8,MPI_MAX,comm_3d,ierr)
         r(np)=tmpout(1)
         Rdelta=-tmpout(2)  
#endif
         !print *,'w,w',winterp,wcontour
         !print *,r(np),Rdelta
         if (Rdelta < contour_eps) exit
         count=count+1
         if (count>5000 .or. r(np)<-9d10) then
            print *,'r(np)=',r(np),np
            print *,'count=',count
            print *,mxcord(1),r(np),cosc(np)
            print *,mxcord(2),r(np),sinc(np)
            call abort("findellipse() count iteration failure")
         endif
      enddo
   enddo
   endif

end subroutine
#endif



#if 0
subroutine findellipse(w,nell,ccord,r)
use params
use mpi
implicit none
real*8 :: w(nx,ny),ccord(2),r(npd)
integer :: nell

!local
real*8 :: mxcord2(2)
real*8 :: Rint(2),wcontour,winterp
real*8 :: Rint_orig(2)
integer :: np,count,ierr
integer :: reset

!
! compute location of ellipse:  ccord(1) + r(np) * cosc(np) 
!                               ccord(2) + r(np) * sinc(np) 
!

   wcontour = wval(nell)*mxw_init
   if (wcontour<mxw) then
   !print *,'looking for ',wcontour,mxw
   do np=1,npd
      ! find Rad so that:  w(Rad cosc(np), Rad sinc(np)) = wcontour
      reset=0
100   continue
      if (r(np)==0) then
         ! first time.  use large min and max values for search  
         Rint(1)=0
         Rint(2)=.9*min(abs(g_xcord(1)-ccord(1)),&
                     abs(g_xcord(o_nx)-ccord(1)),&
                     abs(g_ycord(1)-ccord(2)),&
                     abs(g_ycord(o_ny)-ccord(2)) )
      else
         Rint(1)=r(np)/1.3
         Rint(2)=r(np)*1.3
      endif
      Rint_orig=Rint
      r(np)=(Rint(1)+Rint(2))/2


      count=0
      do 
         call interp4w(w,ccord(1),ccord(2),r(np),cosc(np),sinc(np),winterp)
         if (winterp<=wcontour) Rint(2)=r(np)
         if (winterp>=wcontour)  Rint(1)=r(np)
         if (winterp<-9d10) then
            ! this cpu does not own this point, set data to invalid value
            Rint(1)=-9d20
            Rint(2)=-9d20
         endif
#ifdef USE_MPI
         mxcord2=Rint
         call MPI_allreduce(mxcord2,Rint,2,MPI_REAL8,MPI_MAX,comm_3d,ierr)
#endif
         r(np)=(Rint(1)+Rint(2))/2
         if (Rint(2)-Rint(1) < contour_eps) exit
         count=count+1
         if (count>1000 .or. r(np)<-9d10) call abort("ellipse() count iteration failure")
      enddo


      if (Rint(1)==Rint_orig(1) .or. &
           Rint(2)==Rint_orig(2) ) then
         !print *,'Rint_org: ',Rint_orig
         !print *,'Rint:     ',Rint
         if (reset==1) call abort("ellipse(): bad [rmin,rmax] stopping...")
         call print_message("ellipse(): bad [Rmin,Rmax]. Resetting and trying again...")
         r(np)=0
         reset=1
         goto 100
      endif
   enddo
   endif

end subroutine
#endif




subroutine interp4w(w,xell,yell,r,csc,snc,winterp)
!
!  interpolate to 
!            w(xell + r*cosc,yell + r*sinc) = wc
!
implicit none
real*8 :: r,winterp
real*8 :: w(nx,ny)
real*8 :: xell,yell,csc,snc

!local
real*8 :: x,y,xc,yc
real*8 :: Qint(4)
integer :: jj,igrid,jgrid,jc

x = xell + r*csc
y = yell + r*snc


! interpolate to w(x,y):

! find cpu which owns the grid point x,y
if (xcord(intx1)<=x .and. x<xcord(intx2)+delx .and. &
     ycord(inty1)<=y .and. y<ycord(inty2)+dely ) then

      ! find igrid,jgrid so that point is in box:
      ! igrid-1,igrid,igrid+1,igrid+2   and jgrid-1,jgrid,jgrid+1,jgrid+2
      igrid = intx1 + floor( (x-xcord(intx1))/delx )
      jgrid = inty1 + floor( (y-ycord(inty1))/dely )
      ASSERT("ellipse(): igrid interp error",igrid<=intx2)
      ASSERT("ellipse(): jgrid interp error",jgrid<=inty2)

      ! interpolate trhs
      do jj=1,4
         ! interpolate xcord(igrid-1:igrid+2) to xcord=tracer(i,1)
         ! data  ugrid(igrid-1:igrid+2, jgrid-2+jj,:) 
         xc = 1 + (x-xcord(igrid))/delx
         jc = jgrid-2+jj
         call interp4(w(igrid-1,jc),w(igrid,jc),&
              w(igrid+1,jc),w(igrid+2,jc),&
              xc,Qint(jj))
      enddo
      ! interpolate ycord(jgrid-1:jgrid+2) to ycord=tracer(i,2)
      ! data:  Qint(1:4,j)
      yc = 1 + (y-ycord(jgrid))/dely
      call interp4(Qint(1),Qint(2),Qint(3),Qint(4),yc,winterp)
else
   winterp=-9d20
endif


end subroutine




end module



