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
!
!  read in a DNS .spec file and split it into ASCII gnuplot files
!  one file for r,x,y an z spectrum
!  one file per time snapshot
!
!  to use this program: 
!    set out_base (below) to desired base name of output file
!    set in_name (below) to name of input .spec file
!
!   make splitspec ; ./splitspec
!
!
#include "macros.h"
program splitspec
implicit none
character(len=280) in_name,out_name,out_base,sdata
CPOINTER :: fid
integer ierr,i,n
real*8 :: spec(10000),nr,time


! Output file base name.  
! time and .gp will be appended to the name
out_base = "temp"
in_name = "temp0000.0000.spec"


call copen(in_name,"r",fid,ierr)
if (ierr/=0) then
   print *,'Error opening file: ',in_name
   stop
endif


do
   call cread8e(fid,time,1,ierr)
   print *,'time=',time
   if (ierr/=1) then
      print *,'End of file reached (or error reading input file).  Stopping.'
      exit
   endif

   call cread8(fid,nr,1); n=nr
   print *,'n_r=',n
   call cread8(fid,spec,n)
   
   write(sdata,'(f10.4)') 10000.0000 + time
   out_name = out_base(1:len_trim(out_base)) // sdata(2:10) // "-r.gp"
   open(55,file=out_name)
   do i=1,int(n)
      write(55,'(i4,e14.6)') i-1,spec(i)
   enddo
   close(55)
   
   call cread8(fid,nr,1) ; n=nr
   print *,'n_x=',n
   call cread8(fid,spec,n)
   
   write(sdata,'(f10.4)') 10000.0000 + time
   out_name = out_base(1:len_trim(out_base)) // sdata(2:10) // "-x.gp"
   open(55,file=out_name)
   do i=1,int(n)
      write(55,'(i4,e14.6)') i-1,spec(i)
   enddo
   close(55)
   
   call cread8(fid,nr,1) ; n=nr
   print *,'n_y=',n
   call cread8(fid,spec,n)
   
   write(sdata,'(f10.4)') 10000.0000 + time
   out_name = out_base(1:len_trim(out_base)) // sdata(2:10) // "-y.gp"
   open(55,file=out_name)
   do i=1,int(n)
      write(55,'(i4,e14.6)') i-1,spec(i)
   enddo
   close(55)
   
   call cread8(fid,nr,1) ; n=nr
   print *,'n_z=',n
   call cread8(fid,spec,n)
   
   write(sdata,'(f10.4)') 10000.0000 + time
   out_name = out_base(1:len_trim(out_base)) // sdata(2:10) // "-z.gp"
   open(55,file=out_name)
   do i=1,int(n)
      write(55,'(i4,e14.6)') i-1,spec(i)
   enddo
   close(55)
   
enddo
end program splitspec

