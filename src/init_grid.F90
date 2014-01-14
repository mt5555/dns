!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Copyright 2007.  Los Alamos National Security, LLC. This material was
!produced under U.S. Government contract DE-AC52-06NA25396 for Los
!Alamos National Laboratory (LANL), which is operated by Los Alamos
!National Security, LLC for the U.S. Department of Energy. The
!U.S. Government has rights to use, reproduce, and distribute this
!software.  NEITHER THE GOVERNMENT NOR LOS ALAMOS NATIONAL SECURITY,
!LLC MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY
!FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
!derivative works, such modified software should be clearly marked, so
!as not to confuse it with the version available from LANL.
!
!Additionally, this program is free software; you can redistribute it
!and/or modify it under the terms of the GNU General Public License as
!published by the Free Software Foundation; either version 2 of the
!License, or (at your option) any later version. Accordingly, this
!program is distributed in the hope that it will be useful, but WITHOUT
!ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
!FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
!for more details.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#include "macros.h"
subroutine init_model
use params
use fft_interface
use transpose
implicit none
character(len=80) ::  message
integer :: i,n,k


call params_init()
call transpose_init()
call init_input_file()
call init_grid()     

call mpi_io_init(1) 


if (byteswap_input) call set_byteswap_input(1);


!
! scale alpha now that we know delx
!
if (alpha_value>=1e10) then
   infinite_alpha=1
else if (alpha_value>=1.0) then
   alpha_value=alpha_value*min(delx,dely,delz)
endif
if (infinite_alpha==0) then
   write(message,'(a,f14.8,f10.4)') "NS-Alpha:  alpha, alpha/h :",&
        alpha_value,alpha_value/min(delx,dely,delz)
else
   write(message,'(a,f14.8,f10.4)') "NS-Alpha:  infinite"
endif
call print_message(message)

If (alpha_B>0) then
   write(message,'(a,f14.8)') "Running NS-Bardina-alpha model: alpha_B=",alpha_B
   call print_message(message)
endif


write(message,'(a,i6,a,i6,a,i6)') "Global grid: ",g_nx," x",g_ny," x",g_nz
call print_message(message)
write(message,'(a,i6,a,i6,a,i6)') "Local grid (with padding): ",nx," x",ny," x",nz
call print_message(message)
write(message,'(a,f8.2)') 'kmax/2pi = ',xkmax/2/pi

call fft_interface_init()

!
! are we also tracking passive scalars?
! in NS_UVW model, if n_var > 3, all extra variables are passive
! scalars.
!
! To add passive scalars to the CNS model, we need to update
! cns.F90 to advect them, and here assign all extra variables
! if n_var>5 to be passive scalars
!
if (equations==NS_UVW) then
   ! NS_UVW requires u,v,w (n_var>=3) even for 2d problems
   ! prognostic variables > 3 are passive scalars:
   npassive=max(0,n_var-3)
      
   if (npassive>input_npassive) then
      print *,'WARNING: input file gives data for ',input_npassive,' passive scalars'
      print *,'but dimensions are setup for ',npassive,' passive scalars'
   endif
   if (npassive<input_npassive) then
      if (io_pe==my_pe) then
         print *,'ERROR: input file gives data for ',input_npassive,' passive scalars'
         print *,'but dimensions are setup for ',npassive,' passive scalars'
         call abortdns('stopping...')
      endif
   endif
   npassive=min(npassive,input_npassive)

   if (npassive>0) then
      np1=4
      np2=np1+npassive-1
   endif
endif
write(message,'(a,i6)') "Number of passive scalars to be computed: ",npassive
call print_message(message)

! copy input data into global arrays:
i=0
do n=np1,np2
   i=i+1
   schmidt(n)=input_schmidt(i)
   passive_type(n)=input_passive_type(i)
enddo
if (input_npassive>0) then
deallocate(input_passive_type)
deallocate(input_schmidt)
endif


end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! initialization for fft grid
! domain is from 0..1
! first point is at grid point 0, last point is at gridpoint
! 1-delta_x
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine init_grid
use params
use fft_interface
use mpi
implicit none

!local variables
real*8 :: one=1,xw,maxw
integer i,j,k,l,ierr,i0,i1,j0,j1,k0,k1,ii,jj,iproc_x,iproc_y,splitx
integer im,jm,km
character(len=80) message


if (numerical_method==FOURIER) then
!
! this is needed to gaurentee that the sine and cosine modes
! are on the same processor.  
!
! when fourier computations are done with the reference decomposition:
if (mod(nslabx,2)/=0) then
   call abortdns("nslabx must be even")
endif
if (mod(nslaby,2)/=0) then
   call abortdns("nslaby must be even")
endif
if (mod(nslabz,2)/=0 .and. g_nz>1) then
   call print_message("********************************************************************************")
   call print_message("WARNING:  nslabz must be even if g_nz>1.  Fourier computations done in the ")
   call print_message("          reference decompostion (such as some diagnostics) will be incorrect")
   call print_message("********************************************************************************")
endif

! when fourier computations are done in the z-pencil decomposition
if (mod(nx_2dz,2)/=0) then
   call abortdns("nx_2dz must be even")
endif
if (mod(ny_2dz,2)/=0) then
   call abortdns("ny_2dz must be even")
endif
endif


! periodic FFT case:  for output, we include the point at x=1 (same as x=0)
o_nx=g_nx
if (g_bdy_x1==PERIODIC) o_nx=g_nx+1

o_ny=g_ny
if (g_bdy_y1==PERIODIC) o_ny=g_ny+1

o_nz=g_nz
if (g_bdy_z1==PERIODIC) o_nz=g_nz+1
if (g_nz==1) o_nz=1


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! local grid, bounds over interior (non-boundary) points
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
intx1=nx1
intx2=nx2
inty1=ny1
inty2=ny2
intz1=nz1
intz2=nz2

bx1=nx1
bx2=nx2
by1=ny1
by2=ny2
bz1=nz1
bz2=nz2

if (offset_bdy==1) then
   if (my_x==ncpu_x-1) bx2=nx2+1; o_nx=g_nx+1
   if (my_y==ncpu_y-1) by2=ny2+1; o_ny=g_ny+1
   if (my_z==ncpu_z-1 .and. g_nz>1) then
      bz2=nz2+1; o_nz=g_nz+1
   endif
endif


bdy_x1=g_bdy_x1
bdy_x2=g_bdy_x2
bdy_y1=g_bdy_y1
bdy_y2=g_bdy_y2
bdy_z1=g_bdy_z1
bdy_z2=g_bdy_z2

if (my_x/=0) bdy_x1=INTERNAL
if (my_x/=ncpu_x-1) bdy_x2=INTERNAL
if (my_y/=0) bdy_y1=INTERNAL
if (my_y/=ncpu_y-1) bdy_y2=INTERNAL
if (my_z/=0) bdy_z1=INTERNAL
if (my_z/=ncpu_z-1) bdy_z2=INTERNAL


if REALBOUNDARY(bdy_x1) then
   intx1=bx1+1
endif
if REALBOUNDARY(bdy_x2) then
   intx2=bx2-1
endif
if REALBOUNDARY(bdy_y1) then
   inty1=by1+1
endif
if REALBOUNDARY(bdy_y2) then
   inty2=by2-1
endif
if REALBOUNDARY(bdy_z1) then
   intz1=bz1+1
endif
if REALBOUNDARY(bdy_z2) then
   intz2=bz2-1
endif








!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! global grid data
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
delx = xscale/(o_nx-1)
dely = yscale/(o_ny-1)
if (o_nz==1) then
   delz=zscale
else
   delz = zscale/(o_nz-1)
endif

do i=1,o_nx
   g_xcord(i)=(i-1)*delx	
enddo

do j=1,o_ny
   g_ycord(j)=(j-1)*dely	
enddo
do k=1,o_nz
   g_zcord(k)=(k-1)*delz	
enddo

call fft_get_mcord(g_imcord,g_nx)
call fft_get_mcord(g_jmcord,g_ny)
call fft_get_mcord(g_kmcord,g_nz)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! local grid data  3D decomposition
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

imsine=0              
imsign=0   
imcord=0
do i=bx1,bx2
   l = i-nx1+1 + nslabx*my_x
   xcord(i)=g_xcord(l)
   imcord(i)=g_imcord(l)
   imsine(i)=l-1
   imsign(i)=sign(1,imcord(i))
   if (imcord(i)==0) imsign(i)=0
   if (imcord(i)==g_nx/2) imsign(i)=0
enddo


jmcord=0
jmsign=0
jmsine=0
do j=by1,by2
   l = j-ny1+1 + nslaby*my_y
   ycord(j)=g_ycord(l)
   jmsine(j)=l-1
   jmcord(j)=g_jmcord(l)
   jmsign(j)=sign(1,jmcord(j))
   if (jmcord(j)==0) jmsign(j)=0
   if (jmcord(j)==g_ny/2) jmsign(j)=0
enddo

kmcord=0
kmsign=0
kmsine=0
do k=bz1,bz2
   l = k-nz1+1 + nslabz*my_z
   zcord(k)=g_zcord(l)
   kmcord(k)=g_kmcord(l)
   kmsine(k)=l-1
   kmsign(k)=sign(1,kmcord(k))
   if (kmcord(k)==0) kmsign(k)=0
   if (kmcord(k)==g_nz/2) kmsign(k)=0
enddo


imcord_exp(1:nx)=imcord(1:nx)
jmcord_exp(1:ny)=jmcord(1:ny)
kmcord_exp(1:nz)=kmcord(1:nz)
do i=nx1,nx2,2
   i0=i
   i1=i+1
   if ( imcord(i0)==0 ) imcord_exp(i1)=0
enddo

do j=ny1,ny2,2
   j0=j
   j1=j+1
   if ( jmcord(j0)==0 ) jmcord_exp(j1)=0
enddo

do k=nz1,nz2,2
   k0=k
   k1=k+1
   if ( kmcord(k0)==0 ) kmcord_exp(k1)=0
enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! data for y-pencil decomposition
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do i=1,nx_2dy
   k= -1 + nx1 + i + my_y*nx_2dy    ! nslabx/ncpu_y   
   if (k>nx2) then
      y_imsine(i)=0  ! this is possible if domain decompostion is not perfect
   else
      y_imsine(i)=imsine(k)
   endif
enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! data for z-pencil decomposition
! this code must match what is in transpose_to/from_z()
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
splitx=nslabx/nx_2dz  
iproc_y = my_z/splitx        ! split y direction into this many pieces
iproc_x = mod(my_z,splitx)   ! split x direction into this many pieces

do i=1,nx_2dz
   ii=nx1 + iproc_x*nx_2dz +i -1
   if (ii>nx2 .or. ii<nx1)  then
      call abortdns("Error in z-pencil decomp. setup: bad ii value")
   endif
   z_imcord(i)=imcord(ii)
   z_imsign(i)=imsign(ii)
enddo


do j=1,ny_2dz
   jj = -1 + ny1 + j + iproc_y*ny_2dz      ! ncpy_z*ny_2dz = nslaby
   if (jj>ny2) then
      z_jmcord(j)=0  ! possible, for not-perfect load balance cases
   else
      z_jmcord(j)=jmcord(jj)
   endif
   z_jmsign(j)=sign(1,z_jmcord(j))
   if (z_jmcord(j)==0) z_jmsign(j)=0
   if (z_jmcord(j)==g_ny/2) z_jmsign(j)=0
enddo


z_kmcord=g_kmcord
do k=1,g_nz
   z_kmsign(k)=sign(1,z_kmcord(k))
   if (z_kmcord(k)==0) z_kmsign(k)=0
   if (z_kmcord(k)==g_nz/2) z_kmsign(k)=0
enddo



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! data for x-pencil decomposition
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do k=nz1,nz2
   x_kmcord(k-nz1+1)=kmcord(k)
   x_kmsign(k-nz1+1)=kmsign(k)
enddo

x_imcord=g_imcord
do i=1,g_nx
   x_imsign(i)=sign(1,x_imcord(i))
   if (x_imcord(i)==0) x_imsign(i)=0
   if (x_imcord(i)==g_nx/2) x_imsign(i)=0
enddo


do j=1,ny_2dx
   L = -1 + ny1 + j + my_x*ny_2dx      ! ncpy_z*ny_2dz = nslaby
   if (L>ny2) then
      x_jmcord(j)=0  ! possible, for not-perfect load balance cases
   else
      x_jmcord(j)=jmcord(L)
   endif
   x_jmsign(j)=sign(1,x_jmcord(j))
   if (x_jmcord(j)==0) x_jmsign(j)=0
   if (x_jmcord(j)==g_ny/2) x_jmsign(j)=0
enddo


! compute maximum wave number
maxw=0
do k=nz1,nz2
   km=abs(kmcord(k))
   do j=ny1,ny2
      jm=abs(jmcord(j))
      do i=nx1,nx2
         im=abs(imcord(i))
         xw=im**2 + jm**2 + (km/Lz)**2
         if (dealias_remove(im,jm,km)) xw=0
         maxw=max(maxw,xw)
      enddo
   enddo
enddo
xkmax=2*pi*sqrt(maxw)


end subroutine





subroutine init_input_file()
use params
use transpose
use fft_interface
use mpi
implicit none

!local variables
real*8 :: one=1,xfac,diff_2h,kmode
integer i,j,k,l,ierr
character(len=80) message,carg
integer input_file_type
#ifdef GFORTRAN
integer,intrinsic :: iargc
intrinsic :: getarg
#else
integer,external :: iargc
#endif

if (my_pe==io_pe) then
   !
   ! command line parameters
   ! ./dns -i -r -d rundir  runname  
   !
   !   -d <dirname>  put all output files in dirname/
   !
   !   -i <infile>   take input from infile instead of stdin
   !
   !   -r     use a restart file for the initial conditions
   !   -nors  but dont use a restart file for passive scalars
   !          (usefull for adding scalars to a restart run)
   !
   !   -b   enable bytewsap_input
   !
   !   -nov3   disable ns_vorticity3 code (only in dnsp/ns_xpencil.F90 model) 
   !
   !
   !   -ui  use UDM for input
   !   -uo  use UDM for output
   !
   !   -s   same as -si and -so
   !   -si  input is spectral coefficient file instead of grid point data
   !   -so  output spectral coefficients instead of grid point data
   !   -smax  <n>   output n spectral coefficients
   !
   !   -mio  use mpi-io with a stripe of 64 for grid point I/O
   !
   !   -zi     read input from compressed (NS_UVW) file format
   !   -zo     ouput using compressed (NS_UVW) file format
   !
   !   -cout <text>   type of data to output for the convert.F90 program
   !                   "uvw", "vor", "vorm", "norm2", "passive", "gradu"
   !
   !   -isodir <n>    use at most n directions in isoave calculation
   !   -isodel <n>    use seperation distances at most n grid points.
   !
   !   -i4             set input to 4 bytes
   !   -o4             set output to 4 bytes 
   !
   rundir="./"
   runname="temp"
   i=0
   do 
      i=i+1
      if (i>iargc()) exit
      call getarg(i,message)

      if (message(1:2)=="-b") then
         byteswap_input=.true.
      else if (message(1:5)=="-cout" .and. len_trim(message)==5) then
         i=i+1
         if (i>iargc()) call abortdns("-cout option requires an argument")
         call getarg(i,carg)
         if (carg(1:3)=="uvw") then
            convert_opt=0 
         else if (carg(1:3)=="vor" .and. len_trim(carg)==3) then
            convert_opt=1 
         else if (carg(1:4)=="vorm") then
            convert_opt=2 
         else if (carg(1:5)=="norm2") then
            convert_opt=3 
         else if (carg(1:4)=="4uvw") then
            convert_opt=4 
         else if (carg(1:5)=="norms") then
            convert_opt=5 
         else if (carg(1:7)=="passive") then
            convert_opt=6 
         else if (carg(1:5)=="gradu") then
            convert_opt=7 
         else if (carg(1:15)=="extract_subcube") then
            convert_opt=8 
         else if (carg(1:11)=="spec_window") then
            convert_opt=9 
         else if (carg(1:6)=="iotest") then
            convert_opt=10 
         else if (carg(1:5)=="stats") then
            convert_opt=11 
         else if (carg(1:5)=="uwbar") then
            convert_opt=12 
         else if (carg(1:5)=="trunc") then
            convert_opt=13 
         else if (carg(1:5)=="nlout") then
            convert_opt=14 
         else if (carg(1:4)=="nlin") then
            convert_opt=15 
         else if (carg(1:7)=="dfilter") then
            convert_opt=16 
         else if (carg(1:4)=="dpdf") then
            convert_opt=17 
	 else if (carg(1:5)=="hpass") then
	    convert_opt=18   
	 else if (carg(1:11)=="coarsegrain") then
	    convert_opt=19   
	 else if (carg(1:11)=="dudx") then
	    convert_opt=20
	 else if (carg(1:11)=="ehor") then
	    convert_opt=21
	 else if (carg(1:11)=="potens") then
	    convert_opt=22
	 else if (carg(1:11)=="pe") then
	    convert_opt=23
	 else if (carg(1:11)=="CH_decomp") then
	    convert_opt=24
	 else if (carg(1:6)=="ttrunc") then
            convert_opt=25	          
	 else if (carg(1:6)=="zerocr") then
            convert_opt=26	          
         else
            print *,'cout option: ',carg(1:len_trim(carg))
            call abortdns("-cout unrecognized option")
         endif
      else if (message(1:2)=="-d") then
         i=i+1
         if (i>iargc()) exit
         call getarg(i,rundir)
         j=len_trim(rundir)
         if (j>0) then
            if (rundir(j:j)/="/") then
               rundir(j+1:j+1)="/"
            endif
         endif
      else if (message(1:2)=="-i" .and. len_trim(message)==2) then
         i=i+1
         if (i>iargc()) exit
         call getarg(i,inputfile)
         j=len_trim(inputfile)
      else if (message(1:7)=="-isodir" .and. len_trim(message)==7) then
         i=i+1
         if (i>iargc()) exit
         call getarg(i,carg)
         j=len_trim(carg)
         if (j>0) then
            read(carg,'(i5)') user_specified_isodir
         endif
      else if (message(1:7)=="-isodel" .and. len_trim(message)==7) then
         i=i+1
         if (i>iargc()) exit
         call getarg(i,carg)
         j=len_trim(carg)
         if (j>0) then
            read(carg,'(i5)') user_specified_isodel
         endif
      else if (message(1:4)=="-mio") then
         do_mpi_io=.true.
      else if (message(1:2)=="-r") then
         restart=1
      else if (message(1:5)=="-nors") then
         compute_passive_on_restart=.true.
      else if (message(1:2)=="-s" .and. len_trim(message)==2) then
         r_spec=.true.
         w_spec=.true.
      else if (message(1:3)=="-si" .and. len_trim(message)==3) then
         r_spec=.true.
      else if (message(1:3)=="-so" .and. len_trim(message)==3) then
         w_spec=.true.
      else if (message(1:3)=="-zi" .and. len_trim(message)==3) then
         r_compressed=.true.
      else if (message(1:3)=="-zo" .and. len_trim(message)==3) then
         w_compressed=.true.
      else if (message(1:3)=="-i4" .and. len_trim(message)==3) then
         input_size=4
      else if (message(1:3)=="-o4" .and. len_trim(message)==3) then
         output_size=4
      else if (message(1:5)=="-smax" .and. len_trim(message)==5) then
         i=i+1
         if (i>iargc()) exit
         call getarg(i,carg)
         j=len_trim(carg)
         if (j>0) then
            read(carg,'(i5)') spec_max
         endif
      else if (message(1:3)=="-ui") then
         udm_input=.true.
      else if (message(1:3)=="-uo") then
         udm_output=.true.
      else if (message(1:5)=="-p4wd" .or. message(1:5)=="-p4pg") then
         ! MPICH mpirun adds this "usefull" argument - we have to ignore it
         ! read path:
         i=i+1
         if (i>iargc()) exit
         !call getarg(i,carg)
      else if (message(1:5)=="-nov3") then
         print *,'User disabled ns_vorticity3 code (only if using dnsp model)'
         use_vorticity3=.false.
      else if (message(1:1)/="-") then
         ! this must be the runname
         runname=message(1:len_trim(message))
      else
         print *,'Invalid option (IGNORED):',message
         print *,'./dns [-[t,r,b,s,ui,uo]]  [-i params.inp] [-d rundir] [runname] '
         print *,'-r  restart from rundir/restart.[uvw,h5] '
         print *,'-mio use MPI-IO'
         print *,'more options: see init_grid.F90'
      endif
   enddo
   print *,'Run name:         ',runname(1:len_trim(runname))
   print *,'output directory: ',rundir(1:len_trim(rundir))

   if (len_trim(inputfile)>0) then
      open(5,file=inputfile)
      print *,'reading input from file: ',inputfile	
   endif	
   print *,'Enter input file type: '
   read(5,*) input_file_type
   print *,'input_file_type=',input_file_type	

   if (input_file_type==0) then
      call read_type0()
   else if (input_file_type==2) then
      call read_type2()
   else if (input_file_type==3) then
      call read_type3()
   else if (input_file_type==4) then
      call read_type4()
   else if (input_file_type==5) then
      call read_type5()
   else if (input_file_type==6) then
      call read_type6()
   else if (input_file_type==7) then
      call read_type7()
   else if (input_file_type==8) then
      call read_type8()
   else if (input_file_type==9) then
      call read_type9()
   else if (input_file_type==10) then
      call read_type10()
   else
      ! if you add a new variable and new input type, be sure to add
      ! it to the MPI broadcast!
      call abortdns("bad input file")
   endif


endif


#ifdef USE_MPI
call mpi_bcast(runname,80,MPI_CHARACTER,io_pe,comm_3d ,ierr)
call mpi_bcast(rundir,80,MPI_CHARACTER,io_pe,comm_3d ,ierr)

call mpi_bcast(use_vorticity3,1,MPI_LOGICAL,io_pe,comm_3d ,ierr)
call mpi_bcast(r_spec,1,MPI_LOGICAL,io_pe,comm_3d ,ierr)
call mpi_bcast(w_spec,1,MPI_LOGICAL,io_pe,comm_3d ,ierr)
call mpi_bcast(r_compressed,1,MPI_LOGICAL,io_pe,comm_3d ,ierr)
call mpi_bcast(w_compressed,1,MPI_LOGICAL,io_pe,comm_3d ,ierr)
call mpi_bcast(header_user,1,MPI_INTEGER,io_pe,comm_3d ,ierr)
call mpi_bcast(input_size,1,MPI_INTEGER,io_pe,comm_3d ,ierr)
call mpi_bcast(output_size,1,MPI_INTEGER,io_pe,comm_3d ,ierr)

call mpi_bcast(user_specified_isodel,1,MPI_INTEGER,io_pe,comm_3d ,ierr)
call mpi_bcast(user_specified_isodir,1,MPI_INTEGER,io_pe,comm_3d ,ierr)
call mpi_bcast(spec_max,1,MPI_INTEGER,io_pe,comm_3d ,ierr)
call mpi_bcast(byteswap_input,1,MPI_INTEGER,io_pe,comm_3d ,ierr)
call mpi_bcast(compute_passive_on_restart,1,MPI_INTEGER,io_pe,comm_3d ,ierr)
call mpi_bcast(convert_opt,1,MPI_INTEGER,io_pe,comm_3d ,ierr)
call mpi_bcast(do_mpi_io,1,MPI_INTEGER,io_pe,comm_3d ,ierr)
call mpi_bcast(udm_input,1,MPI_INTEGER,io_pe,comm_3d ,ierr)
call mpi_bcast(udm_output,1,MPI_INTEGER,io_pe,comm_3d ,ierr)
call mpi_bcast(equations,1,MPI_INTEGER,io_pe,comm_3d ,ierr)
call mpi_bcast(mu,1,MPI_REAL8,io_pe,comm_3d ,ierr)
call mpi_bcast(mu_hyper_value,1,MPI_REAL8,io_pe,comm_3d ,ierr)
call mpi_bcast(k_Gt,1,MPI_REAL8,io_pe,comm_3d ,ierr)
call mpi_bcast(mu_hypo_value,1,MPI_REAL8,io_pe,comm_3d ,ierr)
call mpi_bcast(mu_hyper,1,MPI_INTEGER,io_pe,comm_3d ,ierr)
call mpi_bcast(hyper_type,1,MPI_INTEGER,io_pe,comm_3d ,ierr)
call mpi_bcast(mu_hypo,1,MPI_INTEGER,io_pe,comm_3d ,ierr)
call mpi_bcast(dealias,1,MPI_INTEGER,io_pe,comm_3d ,ierr)
call mpi_bcast(use_phaseshift,1,MPI_LOGICAL,io_pe,comm_3d ,ierr)
call mpi_bcast(numerical_method,1,MPI_INTEGER,io_pe,comm_3d ,ierr)
call mpi_bcast(time_final,1,MPI_REAL8,io_pe,comm_3d ,ierr)
call mpi_bcast(cfl_adv,1,MPI_REAL8,io_pe,comm_3d ,ierr)
call mpi_bcast(cfl_vis,1,MPI_REAL8,io_pe,comm_3d ,ierr)
call mpi_bcast(delt_min,1,MPI_REAL8,io_pe,comm_3d ,ierr)
call mpi_bcast(delt_max,1,MPI_REAL8,io_pe,comm_3d ,ierr)
call mpi_bcast(restart_dt,1,MPI_REAL8,io_pe,comm_3d ,ierr)
call mpi_bcast(diag_dt,1,MPI_REAL8,io_pe,comm_3d ,ierr)
call mpi_bcast(model_dt,1,MPI_REAL8,io_pe,comm_3d ,ierr)
call mpi_bcast(screen_dt,1,MPI_REAL8,io_pe,comm_3d ,ierr)
call mpi_bcast(output_dt,1,MPI_REAL8,io_pe,comm_3d ,ierr)
call mpi_bcast(output_vorticity,1,MPI_INTEGER,io_pe,comm_3d ,ierr)
call mpi_bcast(restart,1,MPI_INTEGER,io_pe,comm_3d ,ierr)
call mpi_bcast(init_cond,1,MPI_INTEGER,io_pe,comm_3d ,ierr)
call mpi_bcast(init_cond_subtype,1,MPI_INTEGER,io_pe,comm_3d ,ierr)
call mpi_bcast(init_cond_param1,1,MPI_REAL8,io_pe,comm_3d ,ierr)
call mpi_bcast(init_cond_param2,1,MPI_REAL8,io_pe,comm_3d ,ierr)
call mpi_bcast(input_npassive,1,MPI_INTEGER,io_pe,comm_3d ,ierr)
if (input_npassive>0) then
   if (.not. allocated(input_schmidt)) allocate(input_schmidt(input_npassive))
   if (.not. allocated(input_passive_type)) allocate(input_passive_type(input_npassive))
   call mpi_bcast(input_schmidt,input_npassive,MPI_REAL8,io_pe,comm_3d ,ierr)
   call mpi_bcast(input_passive_type,input_npassive,MPI_INTEGER,io_pe,comm_3d ,ierr)
endif

call mpi_bcast(forcing_type,1,MPI_INTEGER,io_pe,comm_3d ,ierr)
call mpi_bcast(force_theta,1,MPI_LOGICAL,io_pe,comm_3d ,ierr)
call mpi_bcast(forcing_peak_waveno,1,MPI_INTEGER,io_pe,comm_3d ,ierr)
call mpi_bcast(ffval,1,MPI_REAL8,io_pe,comm_3d ,ierr)
call mpi_bcast(fparam1,1,MPI_REAL8,io_pe,comm_3d ,ierr)
call mpi_bcast(g_bdy_x1,1,MPI_INTEGER,io_pe,comm_3d ,ierr)
call mpi_bcast(g_bdy_x2,1,MPI_INTEGER,io_pe,comm_3d ,ierr)
call mpi_bcast(g_bdy_y1,1,MPI_INTEGER,io_pe,comm_3d ,ierr)
call mpi_bcast(g_bdy_y2,1,MPI_INTEGER,io_pe,comm_3d ,ierr)
call mpi_bcast(g_bdy_z1,1,MPI_INTEGER,io_pe,comm_3d ,ierr)
call mpi_bcast(g_bdy_z2,1,MPI_INTEGER,io_pe,comm_3d ,ierr)
call mpi_bcast(diag_struct,1,MPI_INTEGER,io_pe,comm_3d ,ierr)
call mpi_bcast(diag_pdfs,1,MPI_INTEGER,io_pe,comm_3d ,ierr)

call mpi_bcast(fcor,1,MPI_REAL8,io_pe,comm_3d ,ierr)
call mpi_bcast(Lz,1,MPI_REAL8,io_pe,comm_3d ,ierr)
call mpi_bcast(bous,1,MPI_REAL8,io_pe,comm_3d ,ierr)

call mpi_bcast(alpha_value,1,MPI_REAL8,io_pe,comm_3d ,ierr)
call mpi_bcast(alpha_B,1,MPI_REAL8,io_pe,comm_3d ,ierr)
call mpi_bcast(smagorinsky,1,MPI_REAL8,io_pe,comm_3d ,ierr)


call mpi_bcast(ncustom,1,MPI_INTEGER,io_pe,comm_3d ,ierr)
if (.not. allocated(custom)) allocate(custom(ncustom))
if (ncustom>0) then
   call mpi_bcast(custom,ncustom,MPI_REAL8,io_pe,comm_3d,ierr)
endif
call mpi_barrier(comm_3d,ierr)

#endif




! Diffusion of mode k:
! d/dt (.5 u**2)  = - mu 2pi 2pi k k u**2
!
! d/dt (KE) / KE = - 2 mu 2pi 2pi k**2,  
! USE k = the dealiased highest mode, k = g_nx/3
!
! if mu is chosen for a given value of diff_2h,
! 
!    mu = diff_2h / (2 2pi 2pi k**2)
! 


write(message,'(a,e10.4)') 'Diffusion coefficient mu=',mu
call print_message(message)

if (mu>0) then
kmode=1
xfac = 2*2*pi*2*pi*kmode**2
diff_2h =  mu*xfac
write(message,'(a,f5.0,a,f8.2)') 'Diffusion d/dt(KE)/KE on mode k = ',kmode,': ',diff_2h
call print_message(message)

kmode = sqrt( (g_nx**2 + g_ny**2 + g_nz**2)/9.0)
if (ndim==2) kmode = sqrt( (g_nx**2 + g_ny**2 )/4.0)
xfac = 2*2*pi*2*pi*kmode**2
diff_2h =  mu*xfac
write(message,'(a,f5.0,a,f8.2)') 'Diffusion d/dt(KE)/KE on mode k = ',kmode,': ',diff_2h
call print_message(message)
endif


if (mu_hyper_value/=0) then
if (hyper_type == 1) then
   call print_message('Using IMPLICIT hyper diffusion')
endif
write(message,'(a,e10.4,a,i5)') 'Hyper diffusion coefficient mu=',mu_hyper_value, &
'   grad^k, k=',2*mu_hyper
call print_message(message)




!kmode=1
!xfac = 2* (2*pi*2*pi*kmode**2)**mu_hyper
!diff_2h =  mu_hyper_value*xfac
!write(message,'(a,f5.0,a,f8.2)') 'Hyper diffusion d/dt(KE)/KE on mode k = ',kmode,': ',diff_2h
!call print_message(message)

!kmode = sqrt( (g_nx**2 + g_ny**2 + g_nz**2)/9.0)
!if (ndim==2) kmode = sqrt( (g_nx**2 + g_ny**2 )/4.0)
!xfac = 2*(2*pi*2*pi*kmode**2)**mu_hyper
!diff_2h =  mu_hyper_value*xfac
!write(message,'(a,f5.0,a,f8.2)') 'Hyper diffusion d/dt(KE)/KE on mode k = ',kmode,': ',diff_2h
!call print_message(message)

endif

if (mu_hypo/=-1) then
   write(message,'(a,e10.4)') 'Hypo diffusion coefficient mu=',mu_hypo_value
   call print_message(message)

   kmode=1
   xfac = 2* (2*pi*2*pi*kmode**2)**(-mu_hypo)
   diff_2h =  mu_hypo_value*xfac
   write(message,'(a,f5.0,a,f8.2)') 'Hyper diffusion d/dt(KE)/KE on mode k = ',kmode,': ',diff_2h
   call print_message(message)
   
   kmode = sqrt( (g_nx**2 + g_ny**2 + g_nz**2)/9.0)
   if (ndim==2) kmode = sqrt( (g_nx**2 + g_ny**2 )/4.0)
   xfac = 2*(2*pi*2*pi*kmode**2)**(-mu_hypo)
   diff_2h =  mu_hypo_value*xfac
   write(message,'(a,f5.0,a,f8.2)') 'Hyper diffusion d/dt(KE)/KE on mode k = ',kmode,': ',diff_2h
   call print_message(message)
endif

if (fcor/=0) then
   write(message,'(a,e10.4)') 'fcor = ',fcor
   call print_message(message)
endif
if (Lz/=1) then
   write(message,'(a,e10.4)') 'aspect ratio Lz = ',Lz
   call print_message(message)
endif


if (bous/=0) then
   write(message,'(a,e10.4)') 'Boussenesque parameter bous = ',bous
   call print_message(message)
endif


if (cfl_adv<0) then
   cfl_adv=-cfl_adv
   E_kmax_cfl=.true.
endif


end subroutine














subroutine read_type0
use params
implicit none

!local variables
character(len=80) message
character(len=20) sdata
real*8 rvalue
integer i


read(5,'(a)') message
do i=1,80
   if (message(i:i)==' ') then
      runname=message(1:i-1)
      exit
   endif
enddo


read(5,'(a12)') sdata
read(5,*) rvalue
if (sdata=='value') then
   mu=rvalue
else 
   call abortdns("only viscosity type 'value' supported for input format=0")
endif

read(5,'(a12)') sdata
if (sdata=='fft') then
   numerical_method=FOURIER
   dealias=0
else if (sdata=='fft-dealias') then
   numerical_method=FOURIER
   dealias=1
else if (sdata=='fft-sphere') then
   numerical_method=FOURIER
   dealias=2
else if (sdata=='fft-phase') then
   numerical_method=FOURIER
   dealias=2
   use_phaseshift=.true.
else if (sdata=='4th') then
   numerical_method=FOURTH_ORDER
   dealias=0
else
   call abortdns("only 'fft' derivative method supported")
endif


read(5,'(a12)') sdata
if (sdata=='periodic') then
else
   call abortdns("only 'perodic' b.c. supported")
endif

read(5,'(a12)') sdata
if (sdata=='periodic') then
else
   call abortdns("only 'perodic' b.c. supported")
endif

read(5,'(a12)') sdata
if (sdata=='periodic') then
else
   call abortdns("only 'perodic' b.c. supported")
endif

read(5,*) time_final
read(5,*) cfl_adv
read(5,*) cfl_vis
read(5,*) delt_min
read(5,*) delt_max
read(5,*) restart_dt
read(5,*) diag_dt
model_dt=diag_dt
read(5,*) screen_dt
read(5,*) output_dt
read(5,*) ncustom
allocate(custom(ncustom))
do i=1,ncustom
   read(5,*) custom(i)
enddo

init_cond=1           ! KH analytic - default initial condition
init_cond_subtype=0   ! default parameters
forcing_type=0   ! no forcing
diag_struct=0
alpha_value=0


end subroutine








subroutine read_type2
use params
implicit none

!local variables
character(len=80) message
character(len=20) sdata
real*8 rvalue,xfac,kmode
integer i


read(5,'(a12)') sdata
print *,'initial condition: ',sdata
if (sdata=='KH-blob') then
   init_cond=0
else if (sdata=='KH-anal') then
   init_cond=1
else if (sdata=='iso12') then
   init_cond=2
else if (sdata=='sht') then
   init_cond=3
else if (sdata=='vxpair') then
   init_cond=4
else if (sdata=='iso12e') then
   init_cond=5
else 
   print *,'value = >>',sdata,'<<'
   call abortdns("invalid initial condtion specified on line 3 on input file")
endif

read(5,*) init_cond_subtype


read(5,'(a12)') sdata
print *,'forcing: ',sdata
if (sdata=='none') then
   forcing_type=0
else if (sdata=='iso12') then
   forcing_type=1
else if (sdata=='iso12w') then
   forcing_type=2
else 
   call abortdns("invalid forcing type specified on line 4 on input file")
endif


read(5,'(a12)') sdata
read(5,*) rvalue
if (sdata=='value') then
   mu=rvalue
else if (sdata=='kediff') then
   if (ndim==2) then
      kmode=sqrt( (g_nx**2 + g_ny**2 )/4.0)
   else
      kmode=sqrt( (g_nx**2 + g_ny**2 + g_nz**2)/9.0)
   endif
   xfac = 2*2*pi*2*pi*kmode**2
   ! diff_2h =  mu*xfac
   mu = rvalue/xfac
else if (sdata=='hyper') then
   if (ndim==2) then
      kmode=sqrt( (g_nx**2 + g_ny**2 )/4.0)
   else
      kmode=sqrt( (g_nx**2 + g_ny**2 + g_nz**2)/9.0)
   endif
   xfac = 2* (2*pi*kmode)**8
   mu_hyper_value = rvalue/xfac
   mu_hyper = 4
else 
   call abortdns("non supported viscosity type")
endif

read(5,*) alpha_value
read(5,*) diag_struct


read(5,'(a12)') sdata
print *,'method: ',sdata
if (sdata=='fft') then
   numerical_method=FOURIER
   dealias=0
else if (sdata=='fft-dealias') then
   numerical_method=FOURIER
   dealias=1
else if (sdata=='fft-sphere') then
   numerical_method=FOURIER
   dealias=2
else if (sdata=='fft-phase') then
   numerical_method=FOURIER
   dealias=2
   use_phaseshift=.true.
else if (sdata=='4th') then
   numerical_method=FOURTH_ORDER
   dealias=0
else
   print *,'value=',sdata
   call abortdns("only 'fft' derivative method supported")
endif


read(5,'(a12)') sdata
print *,'b.c. x-direction: ',sdata
if (sdata=='periodic') then
   g_bdy_x1=PERIODIC
   g_bdy_x2=PERIODIC
else if (sdata=='in0/onesided') then
   g_bdy_x1=INFLOW0_ONESIDED
   g_bdy_x2=INFLOW0_ONESIDED
else if (sdata=='in0') then
   g_bdy_x1=INFLOW0
   g_bdy_x2=INFLOW0
else
   call abortdns("x boundary condition not supported")
endif

read(5,'(a12)') sdata
print *,'b.c. y-direction: ',sdata
if (sdata=='periodic') then
   g_bdy_y1=PERIODIC
   g_bdy_y2=PERIODIC
else if (sdata=='reflect') then
   g_bdy_y1=REFLECT
   g_bdy_y2=REFLECT
else if (sdata=='custom0') then
   g_bdy_y1=REFLECT_ODD
   g_bdy_y2=INFLOW0_ONESIDED
else if (sdata=='custom1') then
   g_bdy_y1=REFLECT_ODD
   g_bdy_y2=INFLOW0
else
   call abortdns("y boundary condition not supported")
endif

read(5,'(a12)') sdata
print *,'b.c. z-direction: ',sdata
if (sdata=='periodic') then
   g_bdy_z1=PERIODIC
   g_bdy_z2=PERIODIC
else
   call abortdns("only 'perodic' b.c. supported in z-direction")
endif

read(5,*) time_final
read(5,*) cfl_adv
read(5,*) cfl_vis
read(5,*) delt_min
read(5,*) delt_max
read(5,*) restart_dt
read(5,*) diag_dt
model_dt=diag_dt
read(5,*) screen_dt
read(5,*) output_dt
read(5,*) ncustom
allocate(custom(ncustom))
do i=1,ncustom
   read(5,*) custom(i)
enddo


if (numerical_method==FOURIER) then
   ! make sure boundary conditons are periodic:
   if (g_bdy_x1/=PERIODIC .or. g_bdy_x2/=PERIODIC .or. &
       g_bdy_y1/=PERIODIC .or. g_bdy_y2/=PERIODIC .or. &
       g_bdy_z1/=PERIODIC .or. g_bdy_z2/=PERIODIC) then
      call abortdns("FOURIER method requires all boundary conditions be PERIODIC")
   endif
endif


end subroutine








subroutine read_type3
use params
implicit none

!local variables
character(len=80) message
character(len=20) sdata
real*8 rvalue,xfac,kmode
integer i

read(5,'(a12)') sdata
write(*,'(a,a)') 'equations: ',sdata
if (sdata=='ns_uvw') then
   equations=NS_UVW
else if (sdata=='ns_psivor') then
   equations=NS_PSIVOR
else if (sdata=='shallow') then
   equations=SHALLOW
else 
   print *,'value = >>',sdata,'<<'
   call abortdns("invalid equations specified")
endif

read(5,'(a12)') sdata
print *,'initial condition: ',sdata
if (sdata=='KH-blob') then
   init_cond=0
else if (sdata=='KH-anal') then
   init_cond=1
else if (sdata=='iso12') then
   init_cond=2
else if (sdata=='sht') then
   init_cond=3
else if (sdata=='vxpair') then
   init_cond=4
else if (sdata=='iso12e') then
   init_cond=5
else if (sdata=='zero') then
   init_cond=6
else if (sdata=='decay2048') then
   init_cond=7
else if (sdata=='decay2048_e') then
   init_cond=8
else if (sdata=='decay2048_s') then
   init_cond=9
else 
   print *,'value = >>',sdata,'<<'
   call abortdns("invalid initial condtion specified on line 3 on input file")
endif

read(5,*) init_cond_subtype


read(5,'(a12)') sdata
print *,'forcing: ',sdata
if (sdata=='none') then
   forcing_type=0
else if (sdata=='iso12') then
   forcing_type=1
else if (sdata=='iso12w') then
   forcing_type=2
else if (sdata=='iso') then
   forcing_type=3
else if (sdata=='iso23w') then
   forcing_type=4
else if (sdata=='balu') then
   forcing_type=5
else 
   call abortdns("invalid forcing type specified on line 4 on input file")
endif


read(5,'(a12)') sdata
read(5,*) rvalue
if (sdata=='value') then
   mu=rvalue
else if (sdata=='kediff') then
   if (ndim==2) then
      kmode=sqrt( (g_nx**2 + g_ny**2 )/4.0)
   else
      kmode=sqrt( (g_nx**2 + g_ny**2 + g_nz**2)/9.0)
   endif
   xfac = 2*2*pi*2*pi*kmode**2
   ! diff_2h =  mu*xfac
   mu = rvalue/xfac
else if (sdata=='hyper') then
   if (ndim==2) then
      kmode=sqrt( (g_nx**2 + g_ny**2 )/4.0)
   else
      kmode=sqrt( (g_nx**2 + g_ny**2 + g_nz**2)/9.0)
   endif
   xfac = 2* (2*pi*kmode)**8
   mu_hyper_value = rvalue/xfac
   mu_hyper = 4
else 
   call abortdns("non supported viscosity type")
endif

read(5,*) alpha_value
read(5,*) smagorinsky
read(5,*) diag_struct


read(5,'(a12)') sdata
print *,'method: ',sdata
if (sdata=='fft') then
   numerical_method=FOURIER
   dealias=0
else if (sdata=='fft-dealias') then
   numerical_method=FOURIER
   dealias=1
else if (sdata=='fft-sphere') then
   numerical_method=FOURIER
   dealias=2
else if (sdata=='fft-phase') then
   numerical_method=FOURIER
   dealias=2
   use_phaseshift=.true.
else if (sdata=='4th') then
   numerical_method=FOURTH_ORDER
   dealias=0
else
   print *,'value=',sdata
   call abortdns("only 'fft' derivative method supported")
endif


read(5,'(a12)') sdata
print *,'b.c. x-direction: ',sdata
if (sdata=='periodic') then
   g_bdy_x1=PERIODIC
   g_bdy_x2=PERIODIC
else if (sdata=='in0/onesided') then
   g_bdy_x1=INFLOW0_ONESIDED
   g_bdy_x2=INFLOW0_ONESIDED
else if (sdata=='in0') then
   g_bdy_x1=INFLOW0
   g_bdy_x2=INFLOW0
else
   call abortdns("x boundary condition not supported")
endif

read(5,'(a12)') sdata
print *,'b.c. y-direction: ',sdata
if (sdata=='periodic') then
   g_bdy_y1=PERIODIC
   g_bdy_y2=PERIODIC
else if (sdata=='reflect') then
   g_bdy_y1=REFLECT
   g_bdy_y2=REFLECT
else if (sdata=='custom0') then
   g_bdy_y1=REFLECT_ODD
   g_bdy_y2=INFLOW0_ONESIDED
else if (sdata=='custom1') then
   g_bdy_y1=REFLECT_ODD
   g_bdy_y2=INFLOW0
else
   call abortdns("y boundary condition not supported")
endif

read(5,'(a12)') sdata
print *,'b.c. z-direction: ',sdata
if (sdata=='periodic') then
   g_bdy_z1=PERIODIC
   g_bdy_z2=PERIODIC
else
   call abortdns("only 'perodic' b.c. supported in z-direction")
endif

read(5,*) time_final
read(5,*) cfl_adv
read(5,*) cfl_vis
read(5,*) delt_min
read(5,*) delt_max
read(5,*) restart_dt
read(5,*) diag_dt
model_dt=diag_dt
read(5,*) screen_dt
read(5,*) output_dt
read(5,*) ncustom
allocate(custom(ncustom))
do i=1,ncustom
   read(5,*) custom(i)
enddo


if (numerical_method==FOURIER) then
   ! make sure boundary conditons are periodic:
   if (g_bdy_x1/=PERIODIC .or. g_bdy_x2/=PERIODIC .or. &
       g_bdy_y1/=PERIODIC .or. g_bdy_y2/=PERIODIC .or. &
       g_bdy_z1/=PERIODIC .or. g_bdy_z2/=PERIODIC) then
      call abortdns("FOURIER method requires all boundary conditions be PERIODIC")
   endif
endif


end subroutine






subroutine read_type4
use params
implicit none

!local variables
character(len=80) message
character(len=20) sdata
real*8 rvalue,xfac,kmode
integer i

read(5,'(a12)') sdata
write(*,'(a,a)') 'equations: ',sdata
if (sdata=='ns_uvw') then
   equations=NS_UVW
else if (sdata=='ns_psivor') then
   equations=NS_PSIVOR
else if (sdata=='shallow') then
   equations=SHALLOW
else 
   print *,'value = >>',sdata,'<<'
   call abortdns("invalid equations specified")
endif

read(5,'(a12)') sdata
print *,'initial condition: ',sdata
if (sdata=='KH-blob') then
   init_cond=0
else if (sdata=='KH-anal') then
   init_cond=1
else if (sdata=='iso12') then
   init_cond=2
else if (sdata=='sht') then
   init_cond=3
else if (sdata=='vxpair') then
   init_cond=4
else if (sdata=='iso12e') then
   init_cond=5
else if (sdata=='zero') then
   init_cond=6
else if (sdata=='decay2048') then
   init_cond=7
else if (sdata=='decay2048_e') then
   init_cond=8
else if (sdata=='decay2048_s') then
   init_cond=9
else 
   print *,'value = >>',sdata,'<<'
   call abortdns("invalid initial condtion specified on line 3 on input file")
endif

read(5,*) init_cond_subtype


read(5,'(a12)') sdata
print *,'forcing: ',sdata
if (sdata=='none') then
   forcing_type=0
else if (sdata=='iso12') then
   forcing_type=1
else if (sdata=='iso12w') then
   forcing_type=2
else if (sdata=='iso') then
   forcing_type=3
else if (sdata=='iso23w') then
   forcing_type=4
else if (sdata=='balu') then
   forcing_type=5
else if (sdata=='iso12_hel') then
   forcing_type=6
else 
   call abortdns("invalid forcing type specified on line 4 on input file")
endif


!
! viscosity
!
read(5,'(a12)') sdata
read(5,*) rvalue
if (sdata=='value') then
   mu=rvalue
else if (sdata=='kediff') then
   if (ndim==2) then
      kmode=sqrt( (g_nx**2 + g_ny**2 )/4.0)
   else
      kmode=sqrt( (g_nx**2 + g_ny**2 + g_nz**2)/9.0)
   endif
   xfac = 2*2*pi*2*pi*kmode**2
   ! diff_2h =  mu*xfac
   mu = rvalue/xfac
else if (sdata=='hyper') then
   if (ndim==2) then
      kmode=sqrt( (g_nx**2 + g_ny**2 )/4.0)
   else
      kmode=sqrt( (g_nx**2 + g_ny**2 + g_nz**2)/9.0)
   endif
   xfac = 2* (2*pi*kmode)**8
   mu_hyper_value = rvalue/xfac
   mu_hyper = 4
else 
   call abortdns("non supported viscosity type")
endif

!
! hyper viscosity
!
read(5,'(a12)') sdata
read(5,*) rvalue
if (sdata=='none') then
   ! do nothing
else if (sdata=='hyper4') then
   mu_hyper=2
   mu_hyper_value = rvalue
else if (sdata=='hyper8') then
   mu_hyper=4
   mu_hyper_value = rvalue
else 
   call abortdns("non supported hyper viscosity type")
endif

!
! hypo viscosity
!
read(5,'(a12)') sdata
read(5,*) rvalue
if (sdata=='none') then
   ! do nothing
else if (sdata=='hypo2') then
   mu_hypo=1
   mu_hypo_value = rvalue
else 
   call abortdns("non supported hypo viscosity type")
endif




read(5,*) alpha_value
read(5,*) smagorinsky
read(5,*) diag_struct


read(5,'(a12)') sdata
print *,'method: ',sdata
if (sdata=='fft') then
   numerical_method=FOURIER
   dealias=0
else if (sdata=='fft-dealias') then
   numerical_method=FOURIER
   dealias=1
else if (sdata=='fft-sphere') then
   numerical_method=FOURIER
   dealias=2
else if (sdata=='fft-phase') then
   numerical_method=FOURIER
   dealias=2
   use_phaseshift=.true.
else if (sdata=='4th') then
   numerical_method=FOURTH_ORDER
   dealias=0
else
   print *,'value=',sdata
   call abortdns("only 'fft' derivative method supported")
endif


read(5,'(a12)') sdata
print *,'b.c. x-direction: ',sdata
if (sdata=='periodic') then
   g_bdy_x1=PERIODIC
   g_bdy_x2=PERIODIC
else if (sdata=='in0/onesided') then
   g_bdy_x1=INFLOW0_ONESIDED
   g_bdy_x2=INFLOW0_ONESIDED
else if (sdata=='in0') then
   g_bdy_x1=INFLOW0
   g_bdy_x2=INFLOW0
else
   call abortdns("x boundary condition not supported")
endif

read(5,'(a12)') sdata
print *,'b.c. y-direction: ',sdata
if (sdata=='periodic') then
   g_bdy_y1=PERIODIC
   g_bdy_y2=PERIODIC
else if (sdata=='reflect') then
   g_bdy_y1=REFLECT
   g_bdy_y2=REFLECT
else if (sdata=='custom0') then
   g_bdy_y1=REFLECT_ODD
   g_bdy_y2=INFLOW0_ONESIDED
else if (sdata=='custom1') then
   g_bdy_y1=REFLECT_ODD
   g_bdy_y2=INFLOW0
else
   call abortdns("y boundary condition not supported")
endif

read(5,'(a12)') sdata
print *,'b.c. z-direction: ',sdata
if (sdata=='periodic') then
   g_bdy_z1=PERIODIC
   g_bdy_z2=PERIODIC
else
   call abortdns("only 'perodic' b.c. supported in z-direction")
endif

read(5,*) time_final
read(5,*) cfl_adv
read(5,*) cfl_vis
read(5,*) delt_min
read(5,*) delt_max
read(5,*) restart_dt
read(5,*) diag_dt
read(5,*) model_dt
read(5,*) screen_dt
read(5,*) output_dt
read(5,*) ncustom
allocate(custom(ncustom))
do i=1,ncustom
   read(5,*) custom(i)
enddo


if (numerical_method==FOURIER) then
   ! make sure boundary conditons are periodic:
   if (g_bdy_x1/=PERIODIC .or. g_bdy_x2/=PERIODIC .or. &
       g_bdy_y1/=PERIODIC .or. g_bdy_y2/=PERIODIC .or. &
       g_bdy_z1/=PERIODIC .or. g_bdy_z2/=PERIODIC) then
      call abortdns("FOURIER method requires all boundary conditions be PERIODIC")
   endif
endif


end subroutine








subroutine read_type5
use params
implicit none

!local variables
character(len=80) message
character(len=20) sdata
real*8 rvalue,xfac,kmode
integer i

read(5,'(a12)') sdata
write(*,'(a,a)') 'equations: ',sdata
if (sdata=='ns_uvw') then
   equations=NS_UVW
else if (sdata=='ns_psivor') then
   equations=NS_PSIVOR
else if (sdata=='shallow') then
   equations=SHALLOW
else if (sdata(1:3)=='cns') then
   equations=CNS
else 
   print *,'value = >>',sdata,'<<'
   call abortdns("invalid equations specified")
endif

read(5,'(a12)') sdata
print *,'initial condition: ',sdata
if (sdata=='KH-blob') then
   init_cond=0
else if (sdata=='KH-anal') then
   init_cond=1
else if (sdata=='iso12') then
   init_cond=2
else if (sdata=='sht') then
   init_cond=3
else if (sdata=='vxpair') then
   init_cond=4
else if (sdata=='iso12e') then
   init_cond=5
else if (sdata=='zero') then
   init_cond=6
else if (sdata=='decay2048') then
   init_cond=7
else if (sdata=='decay2048_e') then
   init_cond=8
else if (sdata=='decay2048_s') then
   init_cond=9
else if (sdata=='3d_rot') then
   init_cond=10
else 
   print *,'value = >>',sdata,'<<'
   call abortdns("invalid initial condtion specified on line 3 on input file")
endif

read(5,*) init_cond_subtype


read(5,'(a12)') sdata
print *,'forcing: ',sdata
if (sdata=='none') then
   forcing_type=0
else if (sdata=='iso12') then
   forcing_type=1
else if (sdata=='iso12w') then
   forcing_type=2
else if (sdata=='iso') then
   forcing_type=3
else if (sdata=='iso23w') then
   forcing_type=4
else if (sdata=='balu') then
   forcing_type=5
else if (sdata=='iso12_hel') then
   forcing_type=6
else if (sdata=='iso_high_24') then
   forcing_peak_waveno=24
   forcing_type=7
else if (sdata=='iso_high_16') then
   forcing_peak_waveno=16
   forcing_type=7
else if (sdata=='sto_high_24') then
   forcing_peak_waveno=24
   forcing_type=8
else if (sdata=='sto_high_16') then
   forcing_peak_waveno=16
   forcing_type=8
else if (sdata=='sto_high_10') then
   forcing_peak_waveno=10
   forcing_type=8
else 
   call abortdns("invalid forcing type specified on line 4 on input file")
endif


!
! viscosity
!
read(5,'(a12)') sdata
read(5,*) rvalue
if (sdata=='value') then
   mu=rvalue
else if (sdata=='kediff') then
   if (ndim==2) then
      kmode=sqrt( (g_nx**2 + g_ny**2 )/4.0)
   else
      kmode=sqrt( (g_nx**2 + g_ny**2 + g_nz**2)/9.0)
   endif
   xfac = 2*2*pi*2*pi*kmode**2
   ! diff_2h =  mu*xfac
   mu = rvalue/xfac
else if (sdata=='hyper') then
   if (ndim==2) then
      kmode=sqrt( (g_nx**2 + g_ny**2 )/4.0)
   else
      kmode=sqrt( (g_nx**2 + g_ny**2 + g_nz**2)/9.0)
   endif
   xfac = 2* (2*pi*kmode)**8
   mu_hyper_value = rvalue/xfac
   mu_hyper = 4
else 
   call abortdns("non supported viscosity type")
endif

!
! hyper viscosity
!
read(5,'(a12)') sdata
read(5,*) rvalue
if (sdata=='none') then
   ! do nothing
else if (sdata=='hyper4') then
   mu_hyper=2
   mu_hyper_value = rvalue
else if (sdata=='hyper8') then
   mu_hyper=4
   mu_hyper_value = rvalue
else if (sdata=='hyper16') then
   mu_hyper=8
   mu_hyper_value = rvalue
else 
   call abortdns("non supported hyper viscosity type")
endif

!
! hypo viscosity
!
read(5,'(a12)') sdata
read(5,*) rvalue
if (sdata=='none') then
   ! do nothing
else if (sdata=='hypo2') then
   mu_hypo=1
   mu_hypo_value = rvalue
else 
   call abortdns("non supported hypo viscosity type")
endif




read(5,*) alpha_value
read(5,*) smagorinsky
read(5,*) diag_struct
read(5,*) diag_pdfs


read(5,'(a12)') sdata
print *,'method: ',sdata
if (sdata=='fft') then
   numerical_method=FOURIER
   dealias=0
else if (sdata=='fft-dealias') then
   numerical_method=FOURIER
   dealias=1
else if (sdata=='fft-sphere') then
   numerical_method=FOURIER
   dealias=2
else if (sdata=='fft-phase') then
   numerical_method=FOURIER
   dealias=2
   use_phaseshift=.true.
else if (sdata=='4th') then
   numerical_method=FOURTH_ORDER
   dealias=0
else
   print *,'value=',sdata
   call abortdns("only 'fft' derivative method supported")
endif


read(5,'(a12)') sdata
print *,'b.c. x-direction: ',sdata
if (sdata=='periodic') then
   g_bdy_x1=PERIODIC
   g_bdy_x2=PERIODIC
else if (sdata=='in0/onesided') then
   g_bdy_x1=INFLOW0_ONESIDED
   g_bdy_x2=INFLOW0_ONESIDED
else if (sdata=='in0') then
   g_bdy_x1=INFLOW0
   g_bdy_x2=INFLOW0
else
   call abortdns("x boundary condition not supported")
endif

read(5,'(a12)') sdata
print *,'b.c. y-direction: ',sdata
if (sdata=='periodic') then
   g_bdy_y1=PERIODIC
   g_bdy_y2=PERIODIC
else if (sdata=='reflect') then
   g_bdy_y1=REFLECT
   g_bdy_y2=REFLECT
else if (sdata=='custom0') then
   g_bdy_y1=REFLECT_ODD
   g_bdy_y2=INFLOW0_ONESIDED
else if (sdata=='custom1') then
   g_bdy_y1=REFLECT_ODD
   g_bdy_y2=INFLOW0
else
   call abortdns("y boundary condition not supported")
endif

read(5,'(a12)') sdata
print *,'b.c. z-direction: ',sdata
if (sdata=='periodic') then
   g_bdy_z1=PERIODIC
   g_bdy_z2=PERIODIC
else
   call abortdns("only 'perodic' b.c. supported in z-direction")
endif

read(5,*) time_final
read(5,*) cfl_adv
read(5,*) cfl_vis
read(5,*) delt_min
read(5,*) delt_max
read(5,*) restart_dt
read(5,*) diag_dt
read(5,*) model_dt
read(5,*) screen_dt
read(5,*) output_dt
read(5,*) ncustom
allocate(custom(ncustom))
do i=1,ncustom
   read(5,*) custom(i)
enddo


if (numerical_method==FOURIER) then
   ! make sure boundary conditons are periodic:
   if (g_bdy_x1/=PERIODIC .or. g_bdy_x2/=PERIODIC .or. &
       g_bdy_y1/=PERIODIC .or. g_bdy_y2/=PERIODIC .or. &
       g_bdy_z1/=PERIODIC .or. g_bdy_z2/=PERIODIC) then
      call abortdns("FOURIER method requires all boundary conditions be PERIODIC")
   endif
endif


end subroutine




subroutine read_type6
use params
implicit none

!local variables
character(len=80) message
character(len=20) sdata
real*8 rvalue,xfac,kmode
integer i

read(5,'(a12)') sdata
write(*,'(a,a)') 'equations: ',sdata
if (sdata=='ns_uvw') then
   equations=NS_UVW
else if (sdata=='ns_psivor') then
   equations=NS_PSIVOR
else if (sdata=='shallow') then
   equations=SHALLOW
else if (sdata(1:3)=='cns') then
   equations=CNS
else 
   print *,'value = >>',sdata,'<<'
   call abortdns("invalid equations specified")
endif

read(5,'(a12)') sdata
print *,'initial condition: ',sdata
if (sdata=='KH-blob') then
   init_cond=0
else if (sdata=='KH-anal') then
   init_cond=1
else if (sdata=='iso12') then
   init_cond=2
else if (sdata=='sht') then
   init_cond=3
else if (sdata=='vxpair') then
   init_cond=4
else if (sdata=='iso12e') then
   init_cond=5
else if (sdata=='zero') then
   init_cond=6
else if (sdata=='decay2048') then
   init_cond=7
else if (sdata=='decay2048_e') then
   init_cond=8
else if (sdata=='decay2048_s') then
   init_cond=9
else if (sdata=='3d_rot') then
   init_cond=10
else 
   print *,'value = >>',sdata,'<<'
   call abortdns("invalid initial condtion specified on line 3 on input file")
endif

read(5,*) init_cond_subtype


read(5,'(a12)') sdata
print *,'forcing: ',sdata
if (sdata=='none') then
   forcing_type=0
else if (sdata=='iso12') then
   forcing_type=1
else if (sdata=='iso12w') then
   forcing_type=2
else if (sdata=='iso') then
   forcing_type=3
else if (sdata=='iso23w') then
   forcing_type=4
else if (sdata=='balu') then
   forcing_type=5
else if (sdata=='iso12_hel') then
   forcing_type=6
else if (sdata=='iso_high_24') then
   forcing_peak_waveno=24
   forcing_type=7
else if (sdata=='iso_high_16') then
   forcing_peak_waveno=16
   forcing_type=7
else if (sdata=='sto_high_24') then
   forcing_peak_waveno=24
   forcing_type=8
else if (sdata=='sto_high_10') then
   forcing_peak_waveno=10
   forcing_type=8
else 
   call abortdns("invalid forcing type specified on line 4 on input file")
endif


!
! viscosity
!
read(5,'(a12)') sdata
read(5,*) rvalue
if (sdata=='value') then
   mu=rvalue
else if (sdata=='kediff') then
   if (ndim==2) then
      kmode=sqrt( (g_nx**2 + g_ny**2 )/4.0)
   else
      kmode=sqrt( (g_nx**2 + g_ny**2 + g_nz**2)/9.0)
   endif
   xfac = 2*2*pi*2*pi*kmode**2
   ! diff_2h =  mu*xfac
   mu = rvalue/xfac
else if (sdata=='hyper') then
   if (ndim==2) then
      kmode=sqrt( (g_nx**2 + g_ny**2 )/4.0)
   else
      kmode=sqrt( (g_nx**2 + g_ny**2 + g_nz**2)/9.0)
   endif
   xfac = 2* (2*pi*kmode)**8
   mu_hyper_value = rvalue/xfac
   mu_hyper = 4
else 
   call abortdns("non supported viscosity type")
endif

!
! hyper viscosity
!
read(5,'(a12)') sdata
read(5,*) rvalue
if (sdata=='none') then
   ! do nothing
else if (sdata=='hyper4') then
   mu_hyper=2
   mu_hyper_value = rvalue
else if (sdata=='hyper8') then
   mu_hyper=4
   mu_hyper_value = rvalue
else if (sdata=='hyper16') then
   mu_hyper=8
   mu_hyper_value = rvalue
else 
   call abortdns("non supported hyper viscosity type")
endif

!
! hypo viscosity
!
read(5,'(a12)') sdata
read(5,*) rvalue
if (sdata=='none') then
   ! do nothing
else if (sdata=='hypo2') then
   mu_hypo=1
   mu_hypo_value = rvalue
else 
   call abortdns("non supported hypo viscosity type")
endif




read(5,*) alpha_value
read(5,*) smagorinsky
read(5,*) input_npassive
allocate(input_schmidt(input_npassive))
allocate(input_passive_type(input_npassive))
print *,'npassive (from input file): ',input_npassive
do i=1,input_npassive
   read(5,*) input_schmidt(i),input_passive_type(i)
   print *,input_schmidt(i),input_passive_type(i)
enddo



read(5,*) diag_struct
read(5,*) diag_pdfs


read(5,'(a12)') sdata
print *,'method: ',sdata
if (sdata=='fft') then
   numerical_method=FOURIER
   dealias=0
else if (sdata=='fft-dealias') then
   numerical_method=FOURIER
   dealias=1
else if (sdata=='fft-sphere') then
   numerical_method=FOURIER
   dealias=2
else if (sdata=='fft-phase') then
   numerical_method=FOURIER
   dealias=2
   use_phaseshift=.true.
else if (sdata=='4th') then
   numerical_method=FOURTH_ORDER
   dealias=0
else
   print *,'value=',sdata
   call abortdns("only 'fft' derivative method supported")
endif


read(5,'(a12)') sdata
print *,'b.c. x-direction: ',sdata
if (sdata=='periodic') then
   g_bdy_x1=PERIODIC
   g_bdy_x2=PERIODIC
else if (sdata=='in0/onesided') then
   g_bdy_x1=INFLOW0_ONESIDED
   g_bdy_x2=INFLOW0_ONESIDED
else if (sdata=='in0') then
   g_bdy_x1=INFLOW0
   g_bdy_x2=INFLOW0
else
   call abortdns("x boundary condition not supported")
endif

read(5,'(a12)') sdata
print *,'b.c. y-direction: ',sdata
if (sdata=='periodic') then
   g_bdy_y1=PERIODIC
   g_bdy_y2=PERIODIC
else if (sdata=='reflect') then
   g_bdy_y1=REFLECT
   g_bdy_y2=REFLECT
else if (sdata=='custom0') then
   g_bdy_y1=REFLECT_ODD
   g_bdy_y2=INFLOW0_ONESIDED
else if (sdata=='custom1') then
   g_bdy_y1=REFLECT_ODD
   g_bdy_y2=INFLOW0
else
   call abortdns("y boundary condition not supported")
endif

read(5,'(a12)') sdata
print *,'b.c. z-direction: ',sdata
if (sdata=='periodic') then
   g_bdy_z1=PERIODIC
   g_bdy_z2=PERIODIC
else
   call abortdns("only 'perodic' b.c. supported in z-direction")
endif

read(5,*) time_final
read(5,*) cfl_adv
read(5,*) cfl_vis
read(5,*) delt_min
read(5,*) delt_max
read(5,*) restart_dt
read(5,*) diag_dt
read(5,*) model_dt
read(5,*) screen_dt
read(5,*) output_dt
read(5,*) output_vorticity
read(5,*) ncustom
allocate(custom(ncustom))
do i=1,ncustom
   read(5,*) custom(i)
enddo


if (numerical_method==FOURIER) then
   ! make sure boundary conditons are periodic:
   if (g_bdy_x1/=PERIODIC .or. g_bdy_x2/=PERIODIC .or. &
       g_bdy_y1/=PERIODIC .or. g_bdy_y2/=PERIODIC .or. &
       g_bdy_z1/=PERIODIC .or. g_bdy_z2/=PERIODIC) then
      call abortdns("FOURIER method requires all boundary conditions be PERIODIC")
   endif
endif


end subroutine





subroutine read_type7
use params
implicit none

!local variables
character(len=80) message
character(len=20) sdata
real*8 rvalue,xfac,kmode
integer i

read(5,'(i2)') header_user
if (header_user<1 .or. header_user>5) then
   print *,'header_user = ',header_user
   call abortdns("invalid header type specified")
endif

read(5,'(a12)') sdata
write(*,'(a,a)') 'equations: ',sdata
if (sdata=='ns_uvw') then
   equations=NS_UVW
else if (sdata=='ns_psivor') then
   equations=NS_PSIVOR
else if (sdata=='shallow') then
   equations=SHALLOW
else if (sdata(1:3)=='cns') then
   equations=CNS
else 
   print *,'value = >>',sdata,'<<'
   call abortdns("invalid equations specified")
endif

read(5,'(a12)') sdata
print *,'initial condition: ',sdata
if (sdata=='KH-blob') then
   init_cond=0
else if (sdata=='KH-anal') then
   init_cond=1
else if (sdata=='iso12') then
   init_cond=2
else if (sdata=='sht') then
   init_cond=3
else if (sdata=='vxpair') then
   init_cond=4
else if (sdata=='iso12e') then
   init_cond=5
else if (sdata=='zero') then
   init_cond=6
else if (sdata=='decay2048') then
   init_cond=7
else if (sdata=='decay2048_e') then
   init_cond=8
else if (sdata=='decay2048_s') then
   init_cond=9
else if (sdata=='3d_rot') then
   init_cond=10
else 
   print *,'value = >>',sdata,'<<'
   call abortdns("invalid initial condtion specified on line 3 on input file")
endif

read(5,*) init_cond_subtype


read(5,'(a12)') sdata
print *,'forcing: ',sdata
if (sdata=='none') then
   forcing_type=0
else if (sdata=='iso12') then
   forcing_type=1
else if (sdata=='iso12w') then
   forcing_type=2
else if (sdata=='iso') then
   forcing_type=3
else if (sdata=='iso23w') then
   forcing_type=4
else if (sdata=='balu') then
   forcing_type=5
else if (sdata=='iso12_hel') then
   forcing_type=6
else if (sdata=='iso_high_24') then
   forcing_peak_waveno=24
   forcing_type=7
else if (sdata=='iso_high_16') then
   forcing_peak_waveno=16
   forcing_type=7
else if (sdata=='sto_high_24') then
   forcing_peak_waveno=24
   forcing_type=8
else if (sdata=='sto_high_16') then
   forcing_peak_waveno=16
   forcing_type=8
else if (sdata=='sto_high_10') then
   forcing_peak_waveno=10
   forcing_type=8
else if (sdata=='sto_high_8') then
   forcing_peak_waveno=8
   forcing_type=8
else if (sdata=='sto_high_3') then
   forcing_peak_waveno=3
   forcing_type=8
else if (sdata=='sto_high_4') then
   forcing_peak_waveno=4
   forcing_type=8
else 
   call abortdns("invalid forcing type specified on line 4 on input file")
endif


!
! viscosity
!
read(5,'(a12)') sdata
read(5,*) rvalue
if (sdata=='value') then
   mu=rvalue
else if (sdata=='kediff') then
   if (ndim==2) then
      kmode=sqrt( (g_nx**2 + g_ny**2 )/4.0)
   else
      kmode=sqrt( (g_nx**2 + g_ny**2 + g_nz**2)/9.0)
   endif
   xfac = 2*2*pi*2*pi*kmode**2
   ! diff_2h =  mu*xfac
   mu = rvalue/xfac
else if (sdata=='hyper') then
   if (ndim==2) then
      kmode=sqrt( (g_nx**2 + g_ny**2 )/4.0)
   else
      kmode=sqrt( (g_nx**2 + g_ny**2 + g_nz**2)/9.0)
   endif
   xfac = 2* (2*pi*kmode)**8
   mu_hyper_value = rvalue/xfac
   mu_hyper = 4
else 
   call abortdns("non supported viscosity type")
endif

!
! hyper viscosity
!
read(5,'(a12)') sdata
read(5,*) rvalue
if (sdata=='none') then
   ! do nothing
else if (sdata=='hyper0') then
   mu_hyper=0
   mu_hyper_value = rvalue
else if (sdata=='hyper4') then
   mu_hyper=2
   mu_hyper_value = rvalue
else if (sdata=='hyper8') then
   mu_hyper=4
   mu_hyper_value = rvalue
else if (sdata=='hyper16') then
   mu_hyper=8
   mu_hyper_value = rvalue
else 
   call abortdns("non supported hyper viscosity type")
endif

!
! hypo viscosity
!
read(5,'(a12)') sdata
read(5,*) rvalue
if (sdata=='none') then
   ! do nothing
else if (sdata=='hypo2') then
   mu_hypo=1
   mu_hypo_value = rvalue
else 
   call abortdns("non supported hypo viscosity type")
endif

read(5,*) fcor
read(5,*) Lz
read(5,*) bous

read(5,*) alpha_value
read(5,*) smagorinsky
read(5,*) input_npassive
allocate(input_schmidt(input_npassive))
allocate(input_passive_type(input_npassive))
print *,'npassive (from input file): ',input_npassive
do i=1,input_npassive
   read(5,*) input_schmidt(i),input_passive_type(i)
   print *,input_schmidt(i),input_passive_type(i)
enddo



read(5,*) diag_struct
read(5,*) diag_pdfs


read(5,'(a12)') sdata
print *,'method: ',sdata
if (sdata=='fft') then
   numerical_method=FOURIER
   dealias=0
else if (sdata=='fft-dealias') then
   numerical_method=FOURIER
   dealias=1
else if (sdata=='fft-sphere') then
   numerical_method=FOURIER
   dealias=2
else if (sdata=='fft-phase') then
   numerical_method=FOURIER
   dealias=2
   use_phaseshift=.true.
else if (sdata=='4th') then
   numerical_method=FOURTH_ORDER
   dealias=0
else
   print *,'value=',sdata
   call abortdns("only 'fft' derivative method supported")
endif


read(5,'(a12)') sdata
print *,'b.c. x-direction: ',sdata
if (sdata=='periodic') then
   g_bdy_x1=PERIODIC
   g_bdy_x2=PERIODIC
else if (sdata=='in0/onesided') then
   g_bdy_x1=INFLOW0_ONESIDED
   g_bdy_x2=INFLOW0_ONESIDED
else if (sdata=='in0') then
   g_bdy_x1=INFLOW0
   g_bdy_x2=INFLOW0
else
   call abortdns("x boundary condition not supported")
endif

read(5,'(a12)') sdata
print *,'b.c. y-direction: ',sdata
if (sdata=='periodic') then
   g_bdy_y1=PERIODIC
   g_bdy_y2=PERIODIC
else if (sdata=='reflect') then
   g_bdy_y1=REFLECT
   g_bdy_y2=REFLECT
else if (sdata=='custom0') then
   g_bdy_y1=REFLECT_ODD
   g_bdy_y2=INFLOW0_ONESIDED
else if (sdata=='custom1') then
   g_bdy_y1=REFLECT_ODD
   g_bdy_y2=INFLOW0
else
   call abortdns("y boundary condition not supported")
endif

read(5,'(a12)') sdata
print *,'b.c. z-direction: ',sdata
if (sdata=='periodic') then
   g_bdy_z1=PERIODIC
   g_bdy_z2=PERIODIC
else
   call abortdns("only 'perodic' b.c. supported in z-direction")
endif

read(5,*) time_final
read(5,*) cfl_adv
read(5,*) cfl_vis
read(5,*) delt_min
read(5,*) delt_max
read(5,*) restart_dt
read(5,*) diag_dt
read(5,*) model_dt
read(5,*) screen_dt
read(5,*) output_dt
read(5,*) output_vorticity
read(5,*) ncustom
allocate(custom(ncustom))
do i=1,ncustom
   read(5,*) custom(i)
enddo


if (numerical_method==FOURIER) then
   ! make sure boundary conditons are periodic:
   if (g_bdy_x1/=PERIODIC .or. g_bdy_x2/=PERIODIC .or. &
       g_bdy_y1/=PERIODIC .or. g_bdy_y2/=PERIODIC .or. &
       g_bdy_z1/=PERIODIC .or. g_bdy_z2/=PERIODIC) then
      call abortdns("FOURIER method requires all boundary conditions be PERIODIC")
   endif
endif


end subroutine







subroutine read_type8
use params
implicit none

!local variables
character(len=80) message
character(len=20) sdata
real*8 rvalue,xfac,kmode
integer i

read(5,'(i2)') header_user
if (header_user<1 .or. header_user>5) then
   print *,'header_user = ',header_user
   call abortdns("invalid header type specified")
endif

read(5,'(a12)') sdata
write(*,'(a,a)') 'equations: ',sdata
if (sdata=='ns_uvw') then
   equations=NS_UVW
else if (sdata=='ns_psivor') then
   equations=NS_PSIVOR
else if (sdata=='shallow') then
   equations=SHALLOW
else if (sdata(1:3)=='cns') then
   equations=CNS
else 
   print *,'value = >>',sdata,'<<'
   call abortdns("invalid equations specified")
endif

read(5,'(a12)') sdata
print *,'initial condition: ',sdata
if (sdata=='KH-blob') then
   init_cond=0
else if (sdata=='KH-anal') then
   init_cond=1
else if (sdata=='iso12') then
   init_cond=2
else if (sdata=='sht') then
   init_cond=3
else if (sdata=='vxpair') then
   init_cond=4
else if (sdata=='iso12e') then
   init_cond=5
else if (sdata=='zero') then
   init_cond=6
else if (sdata=='decay2048') then
   init_cond=7
else if (sdata=='decay2048_e') then
   init_cond=8
else if (sdata=='decay2048_s') then
   init_cond=9
else if (sdata=='3d_rot') then
   init_cond=10
else if (sdata=='thrmlwnd') then
   init_cond=13
else
   print *,'value = >>',sdata,'<<'
   call abortdns("invalid initial condtion specified on line 3 on input file")
endif

read(5,*) init_cond_subtype


read(5,'(a12)') sdata
print *,'forcing: ',sdata
if (sdata=='none') then
   forcing_type=0
else if (sdata=='iso12') then
   forcing_type=1
else if (sdata=='iso12w') then
   forcing_type=2
else if (sdata=='iso') then
   forcing_type=3
else if (sdata=='iso23w') then
   forcing_type=4
else if (sdata=='balu') then
   forcing_type=5
else if (sdata=='iso12_hel') then
   forcing_type=6
else if (sdata=='iso_high_24') then
   forcing_peak_waveno=24
   forcing_type=7
else if (sdata=='iso_high_16') then
   forcing_peak_waveno=16
   forcing_type=7
else if (sdata=='sto_high_24') then
   forcing_peak_waveno=24
   forcing_type=8
else if (sdata=='sto_high_16') then
   forcing_peak_waveno=16
   forcing_type=8
else if (sdata=='sto_high_10') then
   forcing_peak_waveno=10
   forcing_type=8
else if (sdata=='sto_high_8') then
   forcing_peak_waveno=8
   forcing_type=8
else if (sdata=='sto_high_3') then
   forcing_peak_waveno=3
   forcing_type=8
else if (sdata=='sto_high_4') then
   forcing_peak_waveno=4
   forcing_type=8
else if (sdata=='sto_high_t4') then
   forcing_peak_waveno=4
   forcing_type=8
   force_theta=.true.  ! enable stochastic forcing on THETA (first passive tracer)
else if (sdata=='sto_high_t10') then
   forcing_peak_waveno=10
   forcing_type=8
   force_theta=.true.  ! enable stochastic forcing on THETA (first passive tracer)
else if (sdata=='sto_high_t16') then
   forcing_peak_waveno=16
   forcing_type=8
   force_theta=.true.  ! enable stochastic forcing on THETA (first passive tracer)
else if (sdata=='sto_high_32') then
   forcing_peak_waveno=32
   forcing_type=8
else if (sdata=='sto_high_t32') then
   forcing_peak_waveno=32
   forcing_type=8
   force_theta=.true.  
else if (sdata=='sto_high_64') then
   forcing_peak_waveno=64
   forcing_type=8
else 
   call abortdns("invalid forcing type specified on line 4 on input file")
endif


! read two forcing parameters:
read(5,*) ffval,fparam1

!
! viscosity
!
read(5,'(a12)') sdata
read(5,*) rvalue
if (sdata=='value') then
   mu=rvalue
else if (sdata=='kediff') then
   if (ndim==2) then
      kmode=sqrt( (g_nx**2 + g_ny**2 )/4.0)
   else
      kmode=sqrt( (g_nx**2 + g_ny**2 + g_nz**2)/9.0)
   endif
   xfac = 2*2*pi*2*pi*kmode**2
   ! diff_2h =  mu*xfac
   mu = rvalue/xfac
else if (sdata=='hyper') then
   if (ndim==2) then
      kmode=sqrt( (g_nx**2 + g_ny**2 )/4.0)
   else
      kmode=sqrt( (g_nx**2 + g_ny**2 + g_nz**2)/9.0)
   endif
   xfac = 2* (2*pi*kmode)**8
   mu_hyper_value = rvalue/xfac
   mu_hyper = 4
else 
   call abortdns("non supported viscosity type")
endif

!
! hyper viscosity
!
read(5,'(a12)') sdata
if (sdata=='hyperN') then
   read(5,*) rvalue, mu_hyper
else
   read(5,*) rvalue
endif

if (sdata=='none') then
   ! do nothing
else if (sdata=='hyper0') then
   mu_hyper=0
   mu_hyper_value = rvalue
else if (sdata=='hyper2') then
   mu_hyper=1
   mu_hyper_value = rvalue
else if (sdata=='hyper4') then
   mu_hyper=2
   mu_hyper_value = rvalue
else if (sdata=='hyper6') then
   mu_hyper=3
   mu_hyper_value = rvalue
else if (sdata=='hyper8') then
   mu_hyper=4
   mu_hyper_value = rvalue
else if (sdata=='hyper8_imp') then
   hyper_type=1
   mu_hyper=4
   mu_hyper_value = rvalue
else if (sdata=='hyper10') then
   mu_hyper=5
   mu_hyper_value = rvalue
else if (sdata=='hyper16') then
   mu_hyper=8
   mu_hyper_value = rvalue
else if (sdata=='hyperN') then
   mu_hyper = mu_hyper/2   ! divide by two because we compute (k^2)^mu_hyper
   mu_hyper_value = rvalue
else if (sdata=='hyper2_imp') then
   hyper_type=1
   mu_hyper=1
   mu_hyper_value = rvalue
else if (sdata=='hyper16_imp') then
   hyper_type=1
   mu_hyper=8
   mu_hyper_value = rvalue
else if (sdata=='hyper2_fix') then
   hyper_type=3  ! fixed viscosity coefficient set by user
   mu_hyper=1    ! laplacian^8   grad^16
   mu_hyper_value = rvalue
else if (sdata=='hyper16_fix') then
   hyper_type=3  ! fixed viscosity coefficient set by user
   mu_hyper=8    ! laplacian^8   grad^16
   mu_hyper_value = rvalue
else 
   print *,'value = >>',sdata,'<<'
   call abortdns("non supported hyper viscosity type")
endif

!
! hypo viscosity
!
read(5,'(a12)') sdata
read(5,*) rvalue
if (sdata=='none') then
   ! do nothing
else if (sdata=='hypo2') then
   mu_hypo=1
   mu_hypo_value = rvalue
else if (sdata=='hypo0') then
   mu_hypo=0
   mu_hypo_value = rvalue
else 
   call abortdns("non supported hypo viscosity type")
endif

read(5,*) fcor
read(5,*) Lz
read(5,*) bous

read(5,*) alpha_value
read(5,*) smagorinsky
read(5,*) input_npassive
allocate(input_schmidt(input_npassive))
allocate(input_passive_type(input_npassive))
print *,'npassive (from input file): ',input_npassive
do i=1,input_npassive
   read(5,*) input_schmidt(i),input_passive_type(i)
   print *,input_schmidt(i),input_passive_type(i)
enddo



read(5,*) diag_struct
read(5,*) diag_pdfs


read(5,'(a12)') sdata
print *,'method: ',sdata
if (sdata=='fft') then
   numerical_method=FOURIER
   dealias=0
else if (sdata=='fft-dealias') then
   numerical_method=FOURIER
   dealias=1
else if (sdata=='fft-sphere') then
   numerical_method=FOURIER
   dealias=2
else if (sdata=='fft-23sphere') then
   numerical_method=FOURIER
   dealias=3
else if (sdata=='fft-phase') then
   numerical_method=FOURIER
   dealias=2
   use_phaseshift=.true.
else if (sdata=='4th') then
   numerical_method=FOURTH_ORDER
   dealias=0
else
   print *,'value=',sdata
   call abortdns("only 'fft' derivative method supported")
endif


read(5,'(a12)') sdata
print *,'b.c. x-direction: ',sdata
if (sdata=='periodic') then
   g_bdy_x1=PERIODIC
   g_bdy_x2=PERIODIC
else if (sdata=='in0/onesided') then
   g_bdy_x1=INFLOW0_ONESIDED
   g_bdy_x2=INFLOW0_ONESIDED
else if (sdata=='in0') then
   g_bdy_x1=INFLOW0
   g_bdy_x2=INFLOW0
else
   call abortdns("x boundary condition not supported")
endif

read(5,'(a12)') sdata
print *,'b.c. y-direction: ',sdata
if (sdata=='periodic') then
   g_bdy_y1=PERIODIC
   g_bdy_y2=PERIODIC
else if (sdata=='reflect') then
   g_bdy_y1=REFLECT
   g_bdy_y2=REFLECT
else if (sdata=='custom0') then
   g_bdy_y1=REFLECT_ODD
   g_bdy_y2=INFLOW0_ONESIDED
else if (sdata=='custom1') then
   g_bdy_y1=REFLECT_ODD
   g_bdy_y2=INFLOW0
else
   call abortdns("y boundary condition not supported")
endif

read(5,'(a12)') sdata
print *,'b.c. z-direction: ',sdata
if (sdata=='periodic') then
   g_bdy_z1=PERIODIC
   g_bdy_z2=PERIODIC
else
   call abortdns("only 'perodic' b.c. supported in z-direction")
endif

read(5,*) time_final
read(5,*) cfl_adv
read(5,*) cfl_vis
read(5,*) delt_min
read(5,*) delt_max
read(5,*) restart_dt
read(5,*) diag_dt
read(5,*) model_dt
read(5,*) screen_dt
read(5,*) output_dt
read(5,*) output_vorticity
read(5,*) ncustom
allocate(custom(ncustom))
do i=1,ncustom
   read(5,*) custom(i)
enddo


if (numerical_method==FOURIER) then
   ! make sure boundary conditons are periodic:
   if (g_bdy_x1/=PERIODIC .or. g_bdy_x2/=PERIODIC .or. &
       g_bdy_y1/=PERIODIC .or. g_bdy_y2/=PERIODIC .or. &
       g_bdy_z1/=PERIODIC .or. g_bdy_z2/=PERIODIC) then
      call abortdns("FOURIER method requires all boundary conditions be PERIODIC")
   endif
endif


end subroutine




subroutine read_type9
use params
implicit none

!local variables
character(len=80) message
character(len=20) sdata
real*8 rvalue,xfac,kmode
integer i

read(5,'(i2)') header_user
if (header_user<1 .or. header_user>5) then
   print *,'header_user = ',header_user
   call abortdns("invalid header type specified")
endif

read(5,'(a12)') sdata
write(*,'(a,a)') 'equations: ',sdata
if (sdata=='ns_uvw') then
   equations=NS_UVW
else if (sdata=='ns_psivor') then
   equations=NS_PSIVOR
else if (sdata=='shallow') then
   equations=SHALLOW
else if (sdata(1:3)=='cns') then
   equations=CNS
else 
   print *,'value = >>',sdata,'<<'
   call abortdns("invalid equations specified")
endif

read(5,'(a12)') sdata
print *,'initial condition: ',sdata
if (sdata=='KH-blob') then
   init_cond=0
else if (sdata=='KH-anal') then
   init_cond=1
else if (sdata=='iso12') then
   init_cond=2
else if (sdata=='sht') then
   init_cond=3
else if (sdata=='vxpair') then
   init_cond=4
else if (sdata=='iso12e') then
   init_cond=5
else if (sdata=='zero') then
   init_cond=6
else if (sdata=='decay2048') then
   init_cond=7
else if (sdata=='decay2048_e') then
   init_cond=8
else if (sdata=='decay2048_s') then
   init_cond=9
else if (sdata=='3d_rot') then
   init_cond=10
else if (sdata=='TG') then
   init_cond=11
else if (sdata=='EDQMN') then
   init_cond=12
else 
   print *,'value = >>',sdata,'<<'
   call abortdns("invalid initial condtion specified on line 3 on input file")
endif

read(5,*) init_cond_subtype
read(5,*) init_cond_param1, init_cond_param2


read(5,'(a12)') sdata
print *,'forcing: ',sdata
if (sdata=='none') then
   forcing_type=0
else if (sdata=='iso12') then
   forcing_type=1
else if (sdata=='iso12w') then
   forcing_type=2
else if (sdata=='iso') then
   forcing_type=3
else if (sdata=='iso23w') then
   forcing_type=4
else if (sdata=='balu') then
   forcing_type=5
else if (sdata=='iso12_hel') then
   forcing_type=6
else if (sdata=='fixed1_hel') then
   forcing_type=10
else if (sdata=='iso_high_24') then
   forcing_peak_waveno=24
   forcing_type=7
else if (sdata=='iso_high_16') then
   forcing_peak_waveno=16
   forcing_type=7
else if (sdata=='sto_high_24') then
   forcing_peak_waveno=24
   forcing_type=8
else if (sdata=='sto_high_16') then
   forcing_peak_waveno=16
   forcing_type=8
else if (sdata=='sto_high_10') then
   forcing_peak_waveno=10
   forcing_type=8
else if (sdata=='sto_high_8') then
   forcing_peak_waveno=8
   forcing_type=8
else if (sdata=='sto_high_3') then
   forcing_peak_waveno=3
   forcing_type=8
else if (sdata=='sto_high_4') then
   forcing_peak_waveno=4
   forcing_type=8
else if (sdata=='sto_high_t4') then
   forcing_peak_waveno=4
   forcing_type=8
   force_theta=.true.  ! enable stochastic forcing on THETA (first passive tracer)
else if (sdata=='sto_high_t10') then
   forcing_peak_waveno=10
   forcing_type=8
   force_theta=.true.  ! enable stochastic forcing on THETA (first passive tracer)
else if (sdata=='sto_high_t16') then
   forcing_peak_waveno=16
   forcing_type=8
   force_theta=.true.  ! enable stochastic forcing on THETA (first passive tracer)
else if (sdata=='sto_high_32') then
   forcing_peak_waveno=32
   forcing_type=8
else if (sdata=='sto_high_t32') then
   forcing_peak_waveno=32
   forcing_type=8
   force_theta=.true. 
else if (sdata=='sto_high_64') then
   forcing_peak_waveno=64
   forcing_type=8
else if (sdata=='iso1') then
   forcing_type=9
else 
   call abortdns("invalid forcing type specified on line 4 on input file")
endif


! read two forcing parameters:
read(5,*) ffval,fparam1

!
! viscosity
!
read(5,'(a12)') sdata
read(5,*) rvalue
if (sdata=='value') then
   mu=rvalue
else if (sdata=='kediff') then
   if (ndim==2) then
      kmode=sqrt( (g_nx**2 + g_ny**2 )/4.0)
   else
      kmode=sqrt( (g_nx**2 + g_ny**2 + g_nz**2)/9.0)
   endif
   xfac = 2*2*pi*2*pi*kmode**2
   ! diff_2h =  mu*xfac
   mu = rvalue/xfac
else if (sdata=='hyper') then
   if (ndim==2) then
      kmode=sqrt( (g_nx**2 + g_ny**2 )/4.0)
   else
      kmode=sqrt( (g_nx**2 + g_ny**2 + g_nz**2)/9.0)
   endif
   xfac = 2* (2*pi*kmode)**8
   mu_hyper_value = rvalue/xfac
   mu_hyper = 4
else 
   call abortdns("non supported viscosity type")
endif

!
! hyper viscosity
!
read(5,'(a12)') sdata
if (sdata=='hyperN') then
   read(5,*) rvalue, mu_hyper, k_Gt
else
   read(5,*) rvalue
endif

if (sdata=='none') then
   ! do nothing
else if (sdata=='hyper0') then
   mu_hyper=0
   mu_hyper_value = rvalue
else if (sdata=='hyper2') then
   mu_hyper=1
   mu_hyper_value = rvalue
else if (sdata=='hyper4') then
   mu_hyper=2
   mu_hyper_value = rvalue
else if (sdata=='hyper6') then
   mu_hyper=3
   mu_hyper_value = rvalue
else if (sdata=='hyper8') then
   mu_hyper=4
   mu_hyper_value = rvalue
else if (sdata=='hyper8_imp') then
   hyper_type=1
   mu_hyper=4
   mu_hyper_value = rvalue
else if (sdata=='hyper10') then
   mu_hyper=5
   mu_hyper_value = rvalue
else if (sdata=='hyper16') then
   mu_hyper=8
   mu_hyper_value = rvalue
else if (sdata=='hyperN') then
   mu_hyper = mu_hyper/2   ! divide by two because we compute (k^2)^mu_hyper
!   mu_hyper_value = rvalue*(pi2 * k_Gt)**(-2*mu_hyper)
   mu_hyper_value = rvalue
else if (sdata=='hyper2_imp') then
   hyper_type=1
   mu_hyper=1
   mu_hyper_value = rvalue
else if (sdata=='hyper16_imp') then
   hyper_type=1
   mu_hyper=8
   mu_hyper_value = rvalue
else if (sdata=='hyper2_fix') then
   hyper_type=3  ! fixed viscosity coefficient set by user
   mu_hyper=1    ! laplacian^8   grad^16
   mu_hyper_value = rvalue
else if (sdata=='hyper16_fix') then
   hyper_type=3  ! fixed viscosity coefficient set by user
   mu_hyper=8    ! laplacian^8   grad^16
   mu_hyper_value = rvalue
else 
   call abortdns("non supported hyper viscosity type")
endif

!
! hypo viscosity
!
read(5,'(a12)') sdata
read(5,*) rvalue
if (sdata=='none') then
   ! do nothing
else if (sdata=='hypo2') then
   mu_hypo=1
   mu_hypo_value = rvalue
else if (sdata=='hypo0') then
   mu_hypo=0
   mu_hypo_value = rvalue
else 
   call abortdns("non supported hypo viscosity type")
endif

read(5,*) fcor
read(5,*) Lz
read(5,*) bous

read(5,*) alpha_value
read(5,*) smagorinsky
read(5,*) input_npassive
allocate(input_schmidt(input_npassive))
allocate(input_passive_type(input_npassive))
print *,'npassive (from input file): ',input_npassive
do i=1,input_npassive
   read(5,*) input_schmidt(i),input_passive_type(i)
   print *,input_schmidt(i),input_passive_type(i)
enddo



read(5,*) diag_struct
read(5,*) diag_pdfs


read(5,'(a12)') sdata
print *,'method: ',sdata
if (sdata=='fft') then
   numerical_method=FOURIER
   dealias=0
else if (sdata=='fft-dealias') then
   numerical_method=FOURIER
   dealias=1
else if (sdata=='fft-sphere') then
   numerical_method=FOURIER
   dealias=2
else if (sdata=='fft-23sphere') then
   numerical_method=FOURIER
   dealias=3
else if (sdata=='fft-phase') then
   numerical_method=FOURIER
   dealias=2
   use_phaseshift=.true.
else if (sdata=='4th') then
   numerical_method=FOURTH_ORDER
   dealias=0
else
   print *,'value=',sdata
   call abortdns("only 'fft' derivative method supported")
endif


read(5,'(a12)') sdata
print *,'b.c. x-direction: ',sdata
if (sdata=='periodic') then
   g_bdy_x1=PERIODIC
   g_bdy_x2=PERIODIC
else if (sdata=='in0/onesided') then
   g_bdy_x1=INFLOW0_ONESIDED
   g_bdy_x2=INFLOW0_ONESIDED
else if (sdata=='in0') then
   g_bdy_x1=INFLOW0
   g_bdy_x2=INFLOW0
else
   call abortdns("x boundary condition not supported")
endif

read(5,'(a12)') sdata
print *,'b.c. y-direction: ',sdata
if (sdata=='periodic') then
   g_bdy_y1=PERIODIC
   g_bdy_y2=PERIODIC
else if (sdata=='reflect') then
   g_bdy_y1=REFLECT
   g_bdy_y2=REFLECT
else if (sdata=='custom0') then
   g_bdy_y1=REFLECT_ODD
   g_bdy_y2=INFLOW0_ONESIDED
else if (sdata=='custom1') then
   g_bdy_y1=REFLECT_ODD
   g_bdy_y2=INFLOW0
else
   call abortdns("y boundary condition not supported")
endif

read(5,'(a12)') sdata
print *,'b.c. z-direction: ',sdata
if (sdata=='periodic') then
   g_bdy_z1=PERIODIC
   g_bdy_z2=PERIODIC
else
   call abortdns("only 'perodic' b.c. supported in z-direction")
endif

read(5,*) time_final
read(5,*) cfl_adv
read(5,*) cfl_vis
read(5,*) delt_min
read(5,*) delt_max
read(5,*) restart_dt
read(5,*) diag_dt
read(5,*) model_dt
read(5,*) screen_dt
read(5,*) output_dt
read(5,*) output_vorticity
read(5,*) ncustom
allocate(custom(ncustom))
do i=1,ncustom
   read(5,*) custom(i)
enddo


if (numerical_method==FOURIER) then
   ! make sure boundary conditons are periodic:
   if (g_bdy_x1/=PERIODIC .or. g_bdy_x2/=PERIODIC .or. &
       g_bdy_y1/=PERIODIC .or. g_bdy_y2/=PERIODIC .or. &
       g_bdy_z1/=PERIODIC .or. g_bdy_z2/=PERIODIC) then
      call abortdns("FOURIER method requires all boundary conditions be PERIODIC")
   endif
endif


end subroutine





subroutine read_type10
use params
implicit none

!local variables
character(len=80) message
character(len=20) sdata
real*8 rvalue,xfac,kmode
integer i

read(5,'(i2)') header_user
if (header_user<1 .or. header_user>5) then
   print *,'header_user = ',header_user
   call abortdns("invalid header type specified")
endif

read(5,'(a12)') sdata
write(*,'(a,a)') 'equations: ',sdata
if (sdata=='ns_uvw') then
   equations=NS_UVW
else if (sdata=='ns_psivor') then
   equations=NS_PSIVOR
else if (sdata=='shallow') then
   equations=SHALLOW
else if (sdata(1:3)=='cns') then
   equations=CNS
else 
   print *,'value = >>',sdata,'<<'
   call abortdns("invalid equations specified")
endif

read(5,'(a12)') sdata
print *,'initial condition: ',sdata
if (sdata=='KH-blob') then
   init_cond=0
else if (sdata=='KH-anal') then
   init_cond=1
else if (sdata=='iso12') then
   init_cond=2
else if (sdata=='sht') then
   init_cond=3
else if (sdata=='vxpair') then
   init_cond=4
else if (sdata=='iso12e') then
   init_cond=5
else if (sdata=='zero') then
   init_cond=6
else if (sdata=='decay2048') then
   init_cond=7
else if (sdata=='decay2048_e') then
   init_cond=8
else if (sdata=='decay2048_s') then
   init_cond=9
else if (sdata=='3d_rot') then
   init_cond=10
else if (sdata=='TG') then
   init_cond=11
else if (sdata=='EDQMN') then
   init_cond=12
else 
   print *,'value = >>',sdata,'<<'
   call abortdns("invalid initial condtion specified on line 3 on input file")
endif

read(5,*) init_cond_subtype
read(5,*) init_cond_param1, init_cond_param2


read(5,'(a12)') sdata
print *,'forcing: ',sdata
if (sdata=='none') then
   forcing_type=0
else if (sdata=='iso12') then
   forcing_type=1
else if (sdata=='iso12w') then
   forcing_type=2
else if (sdata=='iso') then
   forcing_type=3
else if (sdata=='iso23w') then
   forcing_type=4
else if (sdata=='balu') then
   forcing_type=5
else if (sdata=='iso12_hel') then
   forcing_type=6
else if (sdata=='fixed1_hel') then
   forcing_type=10
else if (sdata=='iso_high_24') then
   forcing_peak_waveno=24
   forcing_type=7
else if (sdata=='iso_high_16') then
   forcing_peak_waveno=16
   forcing_type=7
else if (sdata=='sto_high_24') then
   forcing_peak_waveno=24
   forcing_type=8
else if (sdata=='sto_high_16') then
   forcing_peak_waveno=16
   forcing_type=8
else if (sdata=='sto_high_10') then
   forcing_peak_waveno=10
   forcing_type=8
else if (sdata=='sto_high_8') then
   forcing_peak_waveno=8
   forcing_type=8
else if (sdata=='sto_high_3') then
   forcing_peak_waveno=3
   forcing_type=8
else if (sdata=='sto_high_4') then
   forcing_peak_waveno=4
   forcing_type=8
else if (sdata=='sto_high_t4') then
   forcing_peak_waveno=4
   forcing_type=8
   force_theta=.true.  ! enable stochastic forcing on THETA (first passive tracer)
else if (sdata=='sto_high_t10') then
   forcing_peak_waveno=10
   forcing_type=8
   force_theta=.true.  ! enable stochastic forcing on THETA (first passive tracer)
else if (sdata=='sto_high_t16') then
   forcing_peak_waveno=16
   forcing_type=8
   force_theta=.true.  ! enable stochastic forcing on THETA (first passive tracer)
else if (sdata=='sto_high_32') then
   forcing_peak_waveno=32
   forcing_type=8
else if (sdata=='sto_high_t32') then
   forcing_peak_waveno=32
   forcing_type=8
   force_theta=.true.  ! enable stochastic forcing on THETA (first passive tracer)
else if (sdata=='sto_high_64') then
   forcing_peak_waveno=64
   forcing_type=8
else if (sdata=='iso1') then
   forcing_type=9
else 
   call abortdns("invalid forcing type specified on line 4 on input file")
endif


! read two forcing parameters:
read(5,*) ffval,fparam1

!
! viscosity
!
read(5,'(a12)') sdata
read(5,*) rvalue
if (sdata=='value') then
   mu=rvalue
else if (sdata=='kediff') then
   if (ndim==2) then
      kmode=sqrt( (g_nx**2 + g_ny**2 )/4.0)
   else
      kmode=sqrt( (g_nx**2 + g_ny**2 + g_nz**2)/9.0)
   endif
   xfac = 2*2*pi*2*pi*kmode**2
   ! diff_2h =  mu*xfac
   mu = rvalue/xfac
else if (sdata=='hyper') then
   if (ndim==2) then
      kmode=sqrt( (g_nx**2 + g_ny**2 )/4.0)
   else
      kmode=sqrt( (g_nx**2 + g_ny**2 + g_nz**2)/9.0)
   endif
   xfac = 2* (2*pi*kmode)**8
   mu_hyper_value = rvalue/xfac
   mu_hyper = 4
else 
   call abortdns("non supported viscosity type")
endif

!
! hyper viscosity
!
read(5,'(a12)') sdata
if (sdata=='hyperN') then
   read(5,*) rvalue, mu_hyper, k_Gt
else
   read(5,*) rvalue
endif

if (sdata=='none') then
   ! do nothing
else if (sdata=='hyper0') then
   mu_hyper=0
   mu_hyper_value = rvalue
else if (sdata=='hyper2') then
   mu_hyper=1
   mu_hyper_value = rvalue
else if (sdata=='hyper4') then
   mu_hyper=2
   mu_hyper_value = rvalue
else if (sdata=='hyper4_fixed') then
   hyper_type=4  ! fixed viscosity coefficient set by user
   mu_hyper=2    ! laplacian^2   grad^4
   mu_hyper_value = rvalue
else if (sdata=='hyper6') then
   mu_hyper=3
   mu_hyper_value = rvalue
else if (sdata=='hyper8') then
   mu_hyper=4
   mu_hyper_value = rvalue
else if (sdata=='hyper8_imp') then
   hyper_type=1
   mu_hyper=4
   mu_hyper_value = rvalue
else if (sdata=='hyper10') then
   mu_hyper=5
   mu_hyper_value = rvalue
else if (sdata=='hyper16') then
   mu_hyper=8
   mu_hyper_value = rvalue
else if (sdata=='hyperN') then
   mu_hyper = mu_hyper/2   ! divide by two because we compute (k^2)^mu_hyper
!   mu_hyper_value = rvalue*(pi2 * k_Gt)**(-2*mu_hyper)
   mu_hyper_value = rvalue
else if (sdata=='hyper2_imp') then
   hyper_type=1
   mu_hyper=1
   mu_hyper_value = rvalue
else if (sdata=='hyper16_imp') then
   hyper_type=1
   mu_hyper=8
   mu_hyper_value = rvalue
else if (sdata=='hyper16_fix') then
   hyper_type=3  ! fixed viscosity coefficient set by user
   mu_hyper=8    ! laplacian^8   grad^16
   mu_hyper_value = rvalue
else 
   call abortdns("non supported hyper viscosity type")
endif

!
! hypo viscosity
!
read(5,'(a12)') sdata
read(5,*) rvalue
if (sdata=='none') then
   ! do nothing
else if (sdata=='hypo2') then
   mu_hypo=1
   mu_hypo_value = rvalue
else if (sdata=='hypo0') then
   mu_hypo=0
   mu_hypo_value = rvalue
else 
   call abortdns("non supported hypo viscosity type")
endif

read(5,*) fcor
read(5,*) Lz
read(5,*) bous

read(5,*) alpha_value
read(5,*) alpha_B
read(5,*) smagorinsky
read(5,*) input_npassive
allocate(input_schmidt(input_npassive))
allocate(input_passive_type(input_npassive))
print *,'npassive (from input file): ',input_npassive
do i=1,input_npassive
   read(5,*) input_schmidt(i),input_passive_type(i)
   print *,input_schmidt(i),input_passive_type(i)
enddo



read(5,*) diag_struct
read(5,*) diag_pdfs


read(5,'(a12)') sdata
print *,'method: ',sdata
if (sdata=='fft') then
   numerical_method=FOURIER
   dealias=0
else if (sdata=='fft-dealias') then
   numerical_method=FOURIER
   dealias=1
else if (sdata=='fft-sphere') then
   numerical_method=FOURIER
   dealias=2
else if (sdata=='fft-23sphere') then
   numerical_method=FOURIER
   dealias=3
else if (sdata=='fft-phase') then
   numerical_method=FOURIER
   dealias=2
   use_phaseshift=.true.
else if (sdata=='4th') then
   numerical_method=FOURTH_ORDER
   dealias=0
else
   print *,'value=',sdata
   call abortdns("only 'fft' derivative method supported")
endif


read(5,'(a12)') sdata
print *,'b.c. x-direction: ',sdata
if (sdata=='periodic') then
   g_bdy_x1=PERIODIC
   g_bdy_x2=PERIODIC
else if (sdata=='in0/onesided') then
   g_bdy_x1=INFLOW0_ONESIDED
   g_bdy_x2=INFLOW0_ONESIDED
else if (sdata=='in0') then
   g_bdy_x1=INFLOW0
   g_bdy_x2=INFLOW0
else
   call abortdns("x boundary condition not supported")
endif

read(5,'(a12)') sdata
print *,'b.c. y-direction: ',sdata
if (sdata=='periodic') then
   g_bdy_y1=PERIODIC
   g_bdy_y2=PERIODIC
else if (sdata=='reflect') then
   g_bdy_y1=REFLECT
   g_bdy_y2=REFLECT
else if (sdata=='custom0') then
   g_bdy_y1=REFLECT_ODD
   g_bdy_y2=INFLOW0_ONESIDED
else if (sdata=='custom1') then
   g_bdy_y1=REFLECT_ODD
   g_bdy_y2=INFLOW0
else
   call abortdns("y boundary condition not supported")
endif

read(5,'(a12)') sdata
print *,'b.c. z-direction: ',sdata
if (sdata=='periodic') then
   g_bdy_z1=PERIODIC
   g_bdy_z2=PERIODIC
else
   call abortdns("only 'perodic' b.c. supported in z-direction")
endif

read(5,*) time_final
read(5,*) cfl_adv
read(5,*) cfl_vis
read(5,*) delt_min
read(5,*) delt_max
read(5,*) restart_dt
read(5,*) diag_dt
read(5,*) model_dt
read(5,*) screen_dt
read(5,*) output_dt
read(5,*) output_vorticity
read(5,*) ncustom
allocate(custom(ncustom))
do i=1,ncustom
   read(5,*) custom(i)
enddo


if (numerical_method==FOURIER) then
   ! make sure boundary conditons are periodic:
   if (g_bdy_x1/=PERIODIC .or. g_bdy_x2/=PERIODIC .or. &
       g_bdy_y1/=PERIODIC .or. g_bdy_y2/=PERIODIC .or. &
       g_bdy_z1/=PERIODIC .or. g_bdy_z2/=PERIODIC) then
      call abortdns("FOURIER method requires all boundary conditions be PERIODIC")
   endif
endif


end subroutine








