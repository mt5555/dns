C
C $Id: fbelimb.f,v 1.1.1.1 2001-03-05 22:42:51 mat Exp $
C
c************************************************************
c
c	Routine Name       - fbelimb
c
c************************************************************
c
c	Computers          - sun, exemplar
c
c	First Version      - September, 1990
c	Revision history   - June, 1992
c                            June, 1997
c
c	Purpose            - Given an LU decomposed matrix in
c                            banded format and a RHS vector,
c                            solve the matrix equation by
c                            successsive forward and backward
c                            elimination
c
c	Usage              - call fbelimb (ib,n,lu,x)
c
c	Arguments
c	   In:     ib      - Bandwidth (assumed odd, centered on diagonal)
c                  n       - Number of equations
c                  lu      - Lower/Upper decomposed matrix (created by ludecb)
c	   In/Out: x       - RHS vector/solution vector
c
c	Required Routines  - None
c
c************************************************************
c
      subroutine fbelimb (ib,n,lu,x)
      implicit none
c
c.....Passed variables
      integer ib                ! Bandwidth (assumed odd, centered on diag)
      integer n                 ! Number of equations/unknowns
      real*8 lu(ib,n)           ! Banded LU decomposed matrix (from ludecb)
      real*8 x(n)               ! RHS vector/solution vector
c
c.....Local variables
      integer ihb               ! Half-bandwidth
      integer id                ! Diagonal index
      integer i                 ! Row index
      integer j                 ! "Full" column index
      integer j1, j2            ! "Full" column index bounds
      integer jhat              ! Banded matrix column index
c
c.....Compute the halfbandwidth and the pivot index
      ihb = (ib-1)/2
      id  = ihb + 1
c
c.....Forward eliminate Lc = b for c (stored temporarily in x)
      do i = 1, n
        j1 = max(1, i-ihb)
        j2 = i-1
        do j = j1, j2
          jhat = j - i + id
          x(i) = x(i) - lu(jhat,i)*x(j)
        end do
        x(i) = x(i) / lu(id,i)
      end do
c
c.....Backward eliminate Ux = c for x (U assumed to have 1's on the diagonal)
      do i = n, 1, -1
        j1 = i+1
        j2 = min(i+ihb, n)
        do j = j1, j2
          jhat = j - i + id
          x(i) = x(i) - lu(jhat,i)*x(j)
        end do
      end do

      return
      end
