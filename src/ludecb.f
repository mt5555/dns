C
C $Id: ludecb.f,v 1.1.1.1 2001-03-05 22:42:51 mat Exp $
C
c************************************************************
c
c	Routine Name       - ludecb
c
c************************************************************
c
c	Computers          - sun, exemplar
c
c	First Version      - September, 1990
c	Latest Revision    - June, 1992
c                            June, 1997
c
c	Purpose            - Given square matrix a in banded
c                            format, compute the lower/upper
c                            decomposition of a (also in banded
c                            format).
c
c	Usage              - call ludecb (ib,n,a,lu)
c
c	Arguments
c	   In:   ib        - Bandwidth (assumed odd, centered about diagonal)
c                n         - Number of equations (rows)
c                a         - Coefficient matrix in banded form
c	   Out:  lu        - Contains both the lower and upper
c                            matrix decomposition of a.  U is
c                            stored in the upper triangular
c                            part of lu and L is stored in the
c                            lower part.  The diagonal belogs to
c                            L and the diagonal of U is assumed
c                            to be all ones.
c
c	Required Routines  - None
c
c************************************************************
c
      subroutine ludecb (ib,n,a,lu)
      implicit none
c
c.....Passed variables
      integer ib                ! Matrix bandwidth
      integer n                 ! Number of equations/unknowns
      real*8 a(ib,n)            ! Banded coefficient matrix
      real*8 lu(ib,n)           ! Banded LU decompostion of matrix a
c
c.....Local variables
      integer ihb               ! Half-bandwidth
      integer id                ! Diagonal index
      integer i                 ! Row index
      integer i1, i2            ! Row index bounds
      integer j                 ! Column index
      integer jihat             ! Banded matrix column index for reduction row
      integer jkhat             ! Banded matrix column index for pivot row
      integer k                 ! Pivot index
      integer k1, k2            ! Pivot index bounds
      integer khat              ! Banded matrix pivot index
c
c.....Compute the half bandwidth and the pivot index
      ihb = (ib-1)/2
      id  = ihb + 1
c
c.....Loop over the columns of lu 
      do j = 1, n
c
c.......Loop over the banded rows of lu
        i1 = max (1, j-ihb)
        i2 = min (j+ihb, n)
        do i = i1, i2
          jihat = j + id - i
c
c.........Compute the lower and upper matrices 
          lu(jihat,i) = a(jihat,i)
          k1 = max(1, i-ihb, j-ihb)
          if (i .lt. j) then
            k2 = i - 1
          else
            k2 = j - 1
          end if
          do k = k1, k2
            jkhat = j + id - k
            khat  = k + id - i
            lu(jihat,i) = lu(jihat,i)-lu(khat,i)*lu(jkhat,k)
          end do
          if (i .lt. j) lu(jihat,i) = lu(jihat,i) / lu(id,i)
        end do
      end do

      return
      end
