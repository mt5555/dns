#include "macros.h"
subroutine time_control(itime,time,Q,q1,q2,q3,work1,work2)
use params
implicit none
real*8 :: Q(nx,ny,nz,n_var)
real*8 :: q1(nx,ny,nz,n_var)
real*8 :: q2(nx,ny,nz,n_var)
real*8 :: q3(nx,ny,nz,n_var)
real*8 :: work1(nx,ny,nz)
real*8 :: work2(nx,ny,nz)
real*8 :: time
integer :: itime

! local variables
integer i,j,k,n
character(len=80) message
character(len=80) fname
real*8 remainder, time_target,mumax, umax,time_next,cfl_used_adv,cfl_used_vis,mx
real*8 divx,divi,tmx1,tmx2,del,delke_tot,lambda
logical,external :: check_time
logical :: doit

real*8,allocatable,save :: ints_save(:,:),maxs_save(:,:)
integer,save :: nscalars=0,nsize=0

call wallclock(tmx1)
time_target=time_final

!
! compute new time step
!
! CFL = delt*umax/delx     delt <= CFL*delx/umax
!
! viscous CFL =  delt*mu/delx^2  delt <= CFL*delx^2/mu
!  
!

umax=maxs(4)
if (g_nz>1) then
   mumax = mu*(1/delx**2 + 1/dely**2 + 1/delz**2)
else
   mumax = mu*(1/delx**2 + 1/dely**2)
endif

delt = cfl_adv/umax                         ! advective CFL
delt = min(delt,cfl_vis/mumax)              ! viscous CFL
delt = max(delt,delt_min)
delt = min(delt,delt_max)


!
!  restart dumps
!
doit=check_time(itime,time,restart_dt,0,0.0,time_next)
time_target=min(time_target,time_next)
if (doit) then
   call multfile_io(time,Q)
endif




!
!  output dumps
!
doit=check_time(itime,time,output_dt,ncustom,custom,time_next)
time_target=min(time_target,time_next)
if (doit) then
   write(message,'(f10.4)') 10000.0000 + time

   fname = rundir(1:len_trim(rundir)) // runname(1:len_trim(runname)) // message(2:10) // ".u"
   call singlefile_io(time,Q(1,1,1,1),fname,work1,work2,0,io_pe)
   fname = rundir(1:len_trim(rundir)) // runname(1:len_trim(runname)) // message(2:10) // ".v"
   call singlefile_io(time,Q(1,1,1,2),fname,work1,work2,0,io_pe)
   fname = rundir(1:len_trim(rundir)) // runname(1:len_trim(runname)) // message(2:10) // ".w"
   call singlefile_io(time,Q(1,1,1,3),fname,work1,work2,0,io_pe)

   call vorticity(q1,Q,work1,work2)
   fname = rundir(1:len_trim(rundir)) // runname(1:len_trim(runname)) // message(2:10) // ".vor"
   call singlefile_io(time,q1(1,1,1,3),fname,work1,work2,0,io_pe)
endif


!
! accumulate scalers into an array, output during 
! diagnositc output
!
nscalars=nscalars+1
if (nscalars > nsize) then
   nsize=nscalars+100
   ! scalar arrays need to be enlarged
   if (allocated(ints_save)) deallocate(ints_save)
   allocate(ints_save(nints,nsize))
   if (allocated(maxs_save)) deallocate(maxs_save)
   allocate(maxs_save(nints,nsize))
endif
ints_save(1:nints,nscalars)=ints(1:nints)
maxs_save(1:nints,nscalars)=maxs(1:nints)



!
! diagnostic output
!
doit=check_time(itime,time,diag_dt,0,0.0,time_next)
time_target=min(time_target,time_next)
if (doit) then
   call output_diags(time,Q,q1,q2,q3,work1,work2,ints_save,maxs_save,nints,nscalars)
   nscalars=0
endif



doit=check_time(itime,time,screen_dt,0,0.0,time_next)
time_target=min(time_target,time_next)
! also output first 5 timesteps, unless screen_dt==0
if (itime<5 .and. screen_dt/=0) doit=.true.



! compute CFL used for next time step.
cfl_used_adv=umax*delt
cfl_used_vis=mumax*delt


!
! display screen output
!
if (doit) then

   write(message,'(a,f9.5,a,i5,a,f9.5)') 'time=',time,'(',itime,')  next output=',time_target
   call print_message(message)	

   write(message,'(a,f9.7,a,f6.3,a,f6.3)') 'for next timestep: delt=',delt,' cfl_adv=',cfl_used_adv,' cfl_vis=',cfl_used_vis
   call print_message(message)	

   write(message,'(a,3f22.15)') 'max: (u,v,w) ',maxs(1),maxs(2),maxs(3)
   call print_message(message)	

   call compute_div(Q,q1,work1,work2,divx,divi)
   write(message,'(3(a,e12.5))') 'max(div)=',divx,'   max(vor)',maxs(5)
   call print_message(message)	

   write(message,'(3(a,e12.5))') '<z-vor>=',ints(4),'   <hel>=',ints(5)
   call print_message(message)	

   lambda=sqrt(  5*(2*ints(1))/ints(2)  )
   write(message,'(3(a,f12.5))') 'R_lambda=',lambda*sqrt(2*ints(1))/mu, &
        '  R=',1/mu
   call print_message(message)	

   write(message,'(a,f13.10,a,f13.4,a,f12.7))') 'Ea: ',&
        ints(1)+.5*alpha_value**2 * ints(2),&
        '   |gradU|: ',ints(2),&    
        '         total d/dt(Ea):',maxs(8)
   call print_message(message)	

   if (alpha_value>0) then
   write(message,'(3(a,f12.7))') 'd/dt(Ea) vis=',&
        -mu*ints(2)-mu*alpha_value**2*ints(10),&
        ' f=',ints(9),'                      tot=',&
        -mu*ints(2)-mu*alpha_value**2*ints(10) + ints(9)
   call print_message(message)	
   endif

   write(message,'(a,f13.10,a,f13.4,a,f12.7)') 'ke: ',ints(1),'  enstropy: ',&
        ints(7),'        total d/dt(ke): ',maxs(9)
   call print_message(message)	
   write(message,'(a,f12.7,a,f12.7,a,f12.7,a,f12.7)') &
     'd/dt(ke) vis=',-mu*ints(2),' f=',ints(3),' alpha=',ints(8),&
     ' total=',(-mu*ints(2)+ints(3)+ints(8))
   call print_message(message)	

   
   call print_message("")
endif


!
! restrict delt so we hit the next time_target
! (do this here so the courant numbers printed above represent the
! values used if the time step is not lowered for output control)
delt = min(delt,time_target-time)



call wallclock(tmx2)
tims(3)=tmx2-tmx1




end subroutine







subroutine multfile_io(time,Q)
use params
implicit none
real*8 :: Q(nx,ny,nz,n_var)
real*8 :: time

! local variables
integer i,j,k,n
real*8 xnx,xny,xnz,xnv
character(len=80) message
character(len=20) tmp
CPOINTER :: fid
integer ierr

n=max(mpidims(1),mpidims(2),mpidims(3))
if (n<10) then
   n=5
else if (n<100) then
   n=4
else if (n<1000) then
   n=3
else if (n<10000) then
   n=2
else 
   call abort("opps, we assumed no more than 10000 cpus along one direction!")
endif
write(tmp,'(i5)') 10000+my_x
message="-" // tmp(n:5)
write(tmp,'(i5)') 10000+my_y
message=message(1:len_trim(message)) // "-" // tmp(n:5)
write(tmp,'(i5)') 10000+my_z
message=message(1:len_trim(message)) // "-" // tmp(n:5) // "-"

write(tmp,'(f10.4)') 10000.0000 + time
message=message(1:len_trim(message)) // tmp(2:10)

message = rundir(1:len_trim(rundir)) // runname(1:len_trim(runname)) // message(1:len_trim(message)) // ".data"

!open(unit=10,file=message,form='binary')
call copen(message,"w",fid,ierr)
if (ierr/=0) then
   write(message,'(a,i5)') "restart_write(): Error opening file errno=",ierr
   call abort(message)
endif



call cwrite8(fid,time,1)
xnv=n_var
xnx=nslabx
xny=nslaby
xnz=nslabz
call cwrite8(fid,xnx,1)
call cwrite8(fid,xny,1)
call cwrite8(fid,xnz,1)
call cwrite8(fid,xnv,1)
do n=1,n_var
do k=nz1,nz2
do j=ny1,ny2
!do i=nx1,nx2
   call cwrite8(fid,Q(nx1,j,k,n),nx2-nx1+1)
!enddo
enddo
enddo
enddo
call cclose(fid)


end subroutine






subroutine output_diags(time,Q,q1,q2,q3,work1,work2,ints_save,maxs_save,nv,nscalars)
use params
use structf
implicit none
real*8 :: Q(nx,ny,nz,n_var)
real*8 :: q1(nx,ny,nz,n_var)
real*8 :: q2(nx,ny,nz,n_var)
real*8 :: q3(nx,ny,nz,n_var)
real*8 :: work1(nx,ny,nz)
real*8 :: work2(nx,ny,nz)
real*8 :: time
integer nv,nscalars
real*8 :: ints_save(nv,nscalars)
real*8 :: maxs_save(nv,nscalars)
integer,parameter :: nints_e=13
real*8 :: ints_e(nints_e)
real*8 :: ints_e2(nints_e)

! local variables
integer i,j,k,n
integer :: iwave,iwave_max,ierr
real*8 spec_x(0:g_nx/2)
real*8 spec_y(0:g_ny/2)
real*8 spec_z(0:g_nz/2)
real*8 x
real*8,allocatable  ::  spectrum(:),spectrum1(:)
character(len=80) :: message
character :: access
CPOINTER fid


! append to output files, unless time=0 create a new file 
access="a"
if (time==0) access="w"

iwave_max=max(g_nx,g_ny,g_nz)
allocate(spectrum(0:iwave_max))
allocate(spectrum1(0:iwave_max))
spectrum=0
spectrum1=0

do i=1,3
   iwave=iwave_max
   call compute_spectrum(Q(:,:,:,i),work1,work2,spectrum1,spec_x,spec_y,spec_z,iwave,io_pe)
   spectrum=spectrum+.5*spectrum1
enddo

write(message,'(a,f10.4)') " KE spectrum",time
call plotASCII(spectrum,iwave,message(1:25))
!call plotASCII(spec_x,g_nx/2,message)
!call plotASCII(spec_y,g_ny/2,message)
!call plotASCII(spec_z,g_nz/2,message)





if (my_pe==io_pe) then
!   write(message,'(f10.4)') 10000.0000 + time
!   message = rundir(1:len_trim(rundir)) // runname(1:len_trim(runname)) // message(2:10) // ".spec"
   message = rundir(1:len_trim(rundir)) // runname(1:len_trim(runname)) // ".spec"
   call copen(message,access,fid,ierr)
   if (ierr/=0) then
      write(message,'(a,i5)') "restart_write(): Error opening file errno=",ierr
      call abort(message)
   endif

   call cwrite8(fid,time,1)
   x=1+iwave; call cwrite8(fid,x,1)
   call cwrite8(fid,spectrum,1+iwave)
   x=1+g_nx/2; call cwrite8(fid,x,1)
   call cwrite8(fid,spec_x,1+g_nx/2)
   x=1+g_ny/2; call cwrite8(fid,x,1)
   call cwrite8(fid,spec_y,1+g_ny/2)
   x=1+g_nz/2; call cwrite8(fid,x,1)
   call cwrite8(fid,spec_z,1+g_nz/2)
   call cclose(fid)
endif
deallocate(spectrum)
deallocate(spectrum1)



!
! output structure functions
!
call compute_all_pdfs(Q,q1,q2,q3,work1,ints_e,nints_e)
if (structf_init==1) then
if (my_pe==io_pe) then
!   write(message,'(f10.4)') 10000.0000 + time
!   message = rundir(1:len_trim(rundir)) // runname(1:len_trim(runname)) // message(2:10) // ".sf"
   message = rundir(1:len_trim(rundir)) // runname(1:len_trim(runname)) // ".sf"
   call copen(message,access,fid,ierr)
   if (ierr/=0) then
      write(message,'(a,i5)') "outputSF(): Error opening file errno=",ierr
      call abort(message)
   endif
endif
call outputSF(time,fid)
if (my_pe==io_pe) call cclose(fid)
endif




if (my_pe==io_pe) then
!   write(message,'(f10.4)') 10000.0000 + time
   message = rundir(1:len_trim(rundir)) // runname(1:len_trim(runname)) // ".scalars"
   call copen(message,access,fid,ierr)
   if (ierr/=0) then
      write(message,'(a,i5)') "diag_output(): Error opening .scalars file errno=",ierr
      call abort(message)
   endif
   x=nv; call cwrite8(fid,x,1)
   x=nscalars; call cwrite8(fid,x,1)
   call cwrite8(fid,mu,1)
   call cwrite8(fid,alpha_value,1)
   call cwrite8(fid,ints_save,nv*nscalars);
   call cwrite8(fid,maxs_save,nv*nscalars);

   x=nints_e; call cwrite8(fid,x,1)
   call cwrite8(fid,time,1)
   call cwrite8(fid,ints_e,nints_e)


   call cclose(fid)
endif



end subroutine









subroutine singlefile_io(time,p,fname,work,work2,read,fpe)
!
! I/O routines where all data goes through a single PE and is
! written to a single file
!
! read=0    write data to file fname
! read=1    read data from file fname
!
! fpe       processor to do the file I/O
!
use params
use transpose
implicit none
integer :: read  ! =1 for read, 0 for write
integer :: fpe
real*8 :: time
real*8 :: p(nx,ny,nz)
real*8 :: work2(nx,ny,nz),work(nx,ny,nz)
character(len=*) :: fname

! local variables
integer i,j,k,n
real*8 xnx,xny,xnz
character(len=80) message
integer n_var_start,ierr
CPOINTER fid


if (my_pe==fpe) then

   if (read==1) then
      call copen(fname,"r",fid,ierr)
      if (ierr/=0) then
         write(message,'(3a,i5)') "singlefile_io(): Read Error file=",fname," error no=",ierr
         call abort(message)
      endif
      call cread8(fid,time,1)
      xnx=o_nx
      xny=o_ny
      xnz=o_nz
      call cread8(fid,xnx,1)
      call cread8(fid,xny,1)
      call cread8(fid,xnz,1)
      if (int(xnx)/=o_nx) call abort("Error: restart file nx <> nx set in params.h");
      if (int(xny)/=o_ny) call abort("Error: restart file ny <> ny set in params.h");
      if (int(xnz)/=o_nz) call abort("Error: restart file nz <> nz set in params.h");
      call cread8(fid,g_xcord(1),o_nx)
      call cread8(fid,g_ycord(1),o_ny)
      call cread8(fid,g_zcord(1),o_nz)
   else
      call copen(fname,"w",fid,ierr)
      if (ierr/=0) then
         write(message,'(3a,i5)') "singlefile_io(): Write Error file=",fname," error no=",ierr
         call abort(message)
      endif
      call cwrite8(fid,time,1)
      xnx=o_nx
      xny=o_ny
      xnz=o_nz
      call cwrite8(fid,xnx,1)
      call cwrite8(fid,xny,1)
      call cwrite8(fid,xnz,1)
      call cwrite8(fid,g_xcord(1),o_nx)
      call cwrite8(fid,g_ycord(1),o_ny)
      call cwrite8(fid,g_zcord(1),o_nz)
   endif
endif

if (read==1) then
   call input1(p,work,work2,fid)
else
   call output1(p,work,work2,fid)
endif
if (my_pe==fpe) call cclose(fid)

end subroutine















logical function check_time(itime,time,dt,ncust,cust,time_next)
use params
implicit none
integer ncust,itime
real*8 :: time,dt,cust(ncust)
real*8 :: time_next  ! output

!local variables
real*8 remainder
integer i
real*8 :: small=1e-7

check_time = .false.
time_next=time_final

! custom times always take precedence:
do i=1,ncust
   if (abs(time-cust(i))<small) then
      check_time=.true.
   else if (time<cust(i)) then
      time_next=min(time_next,cust(i))
      exit 
   endif
enddo

if (dt==0) return


if (time>=time_final-small) check_time=.true.
if (time==0) check_time=.true.


if (dt<0) then
   if (mod(itime,nint(abs(dt)))==0) check_time=.true.
endif

if (dt>0) then
   remainder=time - int((small+time)/dt)*dt
   
   if (remainder < small) then
      check_time=.true.
   endif
   time_next = time-remainder+dt

endif

end function


