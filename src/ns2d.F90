!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  subroutine to compute the RHS of the Navier Stokes equations
!  2D version
!


!  periodic in the j direction
!  wall b.c. in the i direction
!  
!
!  input:  bigq
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine NS(rhs,bigq,t)
implicit none
use params
real*8, intent(out) :: rhs(nvard,imaxd,jmaxd)
real*8, intent(in)  :: bigq(nvard,imaxd,jmaxd)       
real*8, intent(in)  :: t

real*8 ::  p,px,py
real*8 ::  ux,uy,uxx,uyy


do i=imin,imax
do j=jmin,jmax

   call eos(p,bigq(rho_index,i,j))   
   ux =
   uy = 
   px =
   py = 
   uxx = 
   uyy = 

   rhs(u_index,i,j) = &
      -u*ux -v*uy -px/bigq(rho_index,i,j) - viscosity*(uxx+uyy)

   rhs(v_index,i,j) = &
      -u*ux -v*uy -px/bigq(rho_index,i,j) - viscosity*(uxx+uyy)

   rhox = 
   rhoy = 
   rhs(rho_index,i,j) =  

enddo
enddo

! along the boundary, we also have to update the rho and energy equations
do i=imin,imax,(imax-imin)
do j=jmin,jmax,(jmax-jmin)
   rhs(rho_index,i,j) = 
enddo
enddo



end subroutine ns















