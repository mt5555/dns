program stripackd_prb
!
!*******************************************************************************
!
!! STRIPACKD_PRB is a test routine for STRIPACKD.
!
!
!  This driver tests software package STRIPACKD for constructing a 
!  Delaunay triangulation and Voronoi diagram of a set of nodes on 
!  the surface of the unit sphere.
!
!  All STRIPACKD subprograms are tested.
!
!  By default, a triangulation is created from a set of N nodes consisting 
!  of the north pole and N-1 points uniformly distributed around the 
!  60-degree parallel (with constant longitudinal separation).  
!
!  The data is stored as RLAT(I), RLON(I), which are the nodal coordinates 
!  in degrees latitude (-90 to 90) and degrees longitude (-180 to 180).
!
  implicit none
!
  integer, parameter :: nmax = 200
  integer, parameter :: nrow = 9
!
  double precision a
  double precision al
  double precision area
  double precision areas
  double precision ds(nmax)
  double precision elat
  double precision elon
  integer ier
  integer iflag
  logical inside
  integer iwk(2*nmax)
  integer k
  integer ksum
  integer kt
  integer lbtri(6,nmax)
  integer lend(nmax)
  integer lin
  integer list(6*nmax)
  integer listc(6*nmax)
  integer lnew
  integer lp
  integer lpl
  integer, parameter :: lplt = 3
  integer, parameter :: lplv = 4
  integer lptr(6*nmax)
  integer ltri(nrow,2*nmax-4)
  integer lw
  integer n
  integer n0
  integer n1
  integer n2
  integer n3
  integer na
  integer nb
  integer nearnd
  integer nn
  integer nt
  logical numbr
  integer nv
  double precision p(3)
  double precision, parameter :: pltsiz = 7.5D+00
  double precision rc(2*nmax-4)
  double precision rlat(nmax)
  double precision rlon(nmax)
  double precision sc
  character ( len = 80 ) trplot_file_name
  character ( len = 80 ) trplot_title
  double precision v1(3)
  double precision v2(3)
  double precision v3(3)
  double precision vlat
  double precision vlon
  double precision vnrm
  character ( len = 80 ) vrplot_file_name
  character ( len = 80 ) vrplot_title
  double precision x(nmax)
  double precision xc(2*nmax-4)
  double precision y(nmax)
  double precision yc(2*nmax-4)
  double precision z(nmax)
  double precision zc(2*nmax-4)
!
  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'STRIPACKD_PRB'
  write ( *, '(a)' ) '  Tests for STRIPACKD.'
!
!  Generate the default set of nodes as latitudinal and longitudinal
!  coordinates. 
!
  n = 9

  rlat(1) = 90.0D+00
  rlat(2:n) = 60.0D+00

  rlon(1) = 0.0D+00
  do k = 2, n
    rlon(k) = dble ( k - 2 ) * 360.0D+00 / dble ( n - 1 )
  end do

  if ( n < 3 .or. nmax < n ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'STRIPACKD_PRB - Fatal error!'
    write ( *, '(a)' ) '  The value of N is illegal.'
    write ( *, '(a,i6,a)' ) '  3 <= N <= NMAX = ', nmax, ' is required.'
    write ( *, '(a,i6)' ) '  Input N = ', n
    stop
  end if
!
!  Set X and Y to the values of RLON and RLAT, respectively,
!  in radians.  (RLON and RLAT are saved for printing by TRPRNT).
!
  sc = atan ( 1.0D+00 ) / 45.0D+00

  do k = 1, n
    x(k) = sc * rlon(k)
    y(k) = sc * rlat(k)
  end do
!
!  Transform spherical coordinates X and Y to Cartesian
!  coordinates (X,Y,Z) on the unit sphere (X**2 + Y**2 + Z**2 = 1).
!
  call trans ( n, y, x, x, y, z )
!
!  Create the triangulation.
!
  call trmesh ( n, x, y, z, list, lptr, lend, lnew, iwk, iwk(n+1), ds, ier )

  if ( ier == -2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'STRIPACKD_PRB - Warning!'
    write ( *, '(a)' ) '  Error in TRMESH.'
    write ( *, '(a)' ) '  The first three nodes are collinear.'
    stop
  else if ( ier > 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'STRIPACKD_PRB - Fatal error!'
    write ( *, '(a)' ) '  Error in TRMESH.'
    write ( *, '(a)' ) '  Duplicate nodes encountered.'
    stop
  end if
!
!  Print the spherical coordinates and adjacency information.
!
!  IFLAG > 0 indicates that RLON and RLAT only are to be printed.
!
  iflag = 1

  call trprnt ( n, rlon, rlat, z, iflag, list, lptr, lend )
!
!  Test TRLIST and TRLPRT by creating and printing a triangle list.
!
  call trlist ( n, list, lptr, lend, nrow, nt, ltri, ier )

  call trlprt ( n, rlon, rlat, z, iflag, nrow, nt, ltri )
!
!  Test TRPLOT by plotting the portion of the triangulation contained 
!  in the hemisphere centered at E = (ELAT,ELON), where ELAT and ELON
!  are taken to be the center of the range of
!  the nodal latitudes and longitudes.
!
  elat = minval ( rlat(1:n) )
  vlat = maxval ( rlat(1:n) )
  elon = minval ( rlon(1:n) )
  vlon = maxval ( rlon(1:n) )

  elat = ( elat + vlat ) / 2.0D+00
  elon = ( elon + vlon ) / 2.0D+00
  a = 90.0D+00
  numbr = n <= 200

  trplot_title = '(Triangulation created by STRIPACKD_PRB)'

  trplot_file_name = 'stripack_prb_del.eps'

  open ( lplt, file = trplot_file_name )

  call trplot ( lplt, pltsiz, elat, elon, a, n, x, y, z, list, &
    lptr, lend, trplot_title, numbr, ier )

  if ( ier == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TRPLOT created the triangulation plot file: ' // &
      trim ( trplot_file_name )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a,i6)' ) 'TRPLOT returned error code ', ier
  end if
!
!  Test AREAS by computing and printing the area of the
!  convex hull of the nodes (sum of triangle
!  areas) relative to the total surface area (4*Pi).
!
  area = 0.0D+00

  do kt = 1, nt
    n1 = ltri(1,kt)
    n2 = ltri(2,kt)
    n3 = ltri(3,kt)
    v1(1) = x(n1)
    v1(2) = y(n1)
    v1(3) = z(n1)
    v2(1) = x(n2)
    v2(2) = y(n2)
    v2(3) = z(n2)
    v3(1) = x(n3)
    v3(2) = y(n3)
    v3(3) = z(n3)
    area = area + areas ( v1, v2, v3 )
  end do

  area = area / ( 16.0D+00 * atan ( 1.0D+00 ) )

  write ( *, '(a)' ) ' '
  write ( *, '(a,f8.2)' ) '  Relative area of convex hull = ', area
!
!  Test BNODES.  The ordered sequence of boundary nodes is stored in IWK.
!
  call bnodes ( n, list, lptr, lend, iwk, nb, na, nt )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Output from BNODES:'
  write ( *, '(a,i6)' ) '  Number of boundary nodes = ', nb
  write ( *, '(a,i6)' ) '  Number of arcs =           ', na
  write ( *, '(a,i6)' ) '  Number of triangles =      ', nt
!
!  Test GETNP by ordering the nodes on distance from N0 and verifying 
!  the ordering.  
!
!  The sequence of nodal indexes is stored in IWK, and the values of an
!  increasing function (the negative cosine) of angular distance is
!  stored in DS.
!
  n0 = n / 2
  iwk(1) = n0
  ds(1) = -1.0D+00
  ksum = n0

  do k = 2, n

    call getnp ( x, y, z, list, lptr, lend, k, iwk, ds(k), ier )

    if ( ier /= 0  .or.  ds(k) < ds(k-1) ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'STRIPACKD_PRB - Fatal error!'
      write ( *, '(a)' ) '  Error in GETNP.'
      stop
    end if

    ksum = ksum + iwk(k)

  end do
!
!  Test for all nodal indexes included in IWK.
!
  if ( ksum /= ( n * ( n + 1 ) ) / 2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'STRIPACKD_PRB - Fatal error!'
    write ( *, '(a)' ) '  Error in GETNP.'
    stop
  end if
!
!  Test NEARND by verifying that the nearest node to K is
!  node K for K = 1 to N.
!
  do k = 1, n

    p(1) = x(k)
    p(2) = y(k)
    p(3) = z(k)

    n0 = nearnd ( p, 1, n, x, y, z, list, lptr, lend, al )

    if ( n0 /= k .or. al > 0.001D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'STRIPACKD_PRB - Fatal error!'
      write ( *, '(a)' ) '  Error in NEARND.'
      stop
    end if

  end do
!
!  Test DELARC by removing a boundary arc if possible.
!  The last two nodes define a boundary arc
!  in the default data set.
!
  n1 = n-1
  n2 = n
  call delarc ( n, n1, n2, list, lptr, lend, lnew, ier )

  if ( ier == 1  .or.  ier == 4 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'STRIPACKD_PRB - Warning!'
    write ( *, '(a,i6)' ) '  DELARC returned error code ', ier
    stop
  end if

  if ( ier /= 0 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Subroutine DELARC was not tested.'
    write ( *, '(a,i6,a,i6,a)' ) '  Nodes ', n1, ' and ', n2, &
      ' do not form a removable boundary arc.'
  else

    call trmesh ( n, x, y, z, list, lptr, lend, lnew, iwk, iwk(n+1), ds, &
      ier )

  end if
!
!  Test CRLIST, VRPLOT, and SCOORD by constructing and
!  plotting the Voronoi diagram, and printing
!  the Voronoi region boundary (ordered
!  sequence of Voronoi vertices) associated with N0.
!
!  Note that the triangulation data structure
!  is altered if NB > 0.
!
  call crlist ( n, nmax, x, y, z, list, lend, lptr, lnew, &
    lbtri, listc, nb, xc, yc, zc, rc, ier )

  if ( ier /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'STRIPACKD_PRB - Warning!'
    write ( *, '(a,i6)' ) '  CRLIST returned error code ', ier
    stop
  end if
!
!  Use the same parameter values that were used for the
!  triangulation plot (except the output unit and title).
!
  nt = 2 * n - 4

  vrplot_file_name = 'stripack_prb_vor.eps'

  vrplot_title = '(Voronoi diagram created by STRIPACKD_PRB)'

  open ( lplv, file = vrplot_file_name )

  call vrplot ( lplv, pltsiz, elat, elon, a, n, x, y, z, nt, listc, &
    lptr, lend, xc, yc, zc, vrplot_title, numbr, ier )

  if ( ier == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'VRPLOT created the Voronoi plot file: ' // &
      trim ( vrplot_file_name )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'STRIPACKD_PRB - Warning!'
    write ( *, '(a,i6)' ) '  VRPLOT returned error code ', ier
  end if

  n0 = 1

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Voronoi region for node ', n0
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     Triangle     Latitude     Longitude' // &
    '     Circumradius'
  write ( *, '(a)' ) ' '
!
!  Initialize for loop on Voronoi vertices (triangle circumcenters).  
!  The number of vertices is accumulated in NV, and the vertex indexes
!  are stored in IWK.  The vertices are converted to latitude and longitude 
!  in degrees for printing.
!
  nv = 0
  lpl = lend(n0)
  lp = lpl

  do

    lp = lptr(lp)
    kt = listc(lp)
    nv = nv + 1
    iwk(nv) = kt
    call scoord ( xc(kt), yc(kt), zc(kt), vlat, vlon, vnrm )
    vlat = vlat / sc
    vlon = vlon / sc
    write ( *, '(i13,f13.6,f14.6,f17.6)' ) kt, vlat, vlon, rc(kt)

    if ( lp == lpl ) then
      exit
    end if

  end do
!
!  Test INSIDE by checking for node N0 inside its Voronoi region.
!
  p(1) = x(n0)
  p(2) = y(n0)
  p(3) = z(n0)

  if ( .not. inside ( p, 2*nmax-4, xc, yc, zc, nv, iwk, ier ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'STRIPACKD_PRB - Warning!'
    write ( *, '(a)' ) '  Error in INSIDE.'
    write ( *, '(a)' ) '  A node is not contained in its Voronoi region.'
    write ( *, '(a,i6)' ) '  Node index = ', n0
  end if

  if ( ier /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'STRIPACKD_PRB - Fatal error!'
    write ( *, '(a)' ) '  Error in INSIDE.'
    write ( *, '(a,i6)' ) '  IER = ', ier
    stop
  end if
!
!  Recreate the triangulation and test the error flag.
!
  call trmesh ( n, x, y, z, list, lptr, lend, lnew, iwk, iwk(n+1), ds, ier )

  if ( ier == -2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'STRIPACKD_PRB - Warning!'
    write ( *, '(a)' ) '  Error in TRMESH.'
    write ( *, '(a)' ) '  The first three nodes are collinear.'
    stop
  else if ( ier > 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'STRIPACKD_PRB - Fatal error!'
    write ( *, '(a)' ) '  Error in TRMESH.'
    write ( *, '(a)' ) '  Duplicate nodes encountered.'
    stop
  end if
!
!  Test EDGE by forcing an edge between nodes N1=1 and N2=N. 
!
  n1 = 1
  n2 = n

  call edge ( n1, n2, x, y, z, nmax, iwk, list, lptr, lend, ier )

  if ( ier /= 0 .and. ier /= 5 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'STRIPACKD_PRB - Fatal error!'
    write ( *, '(a)' ) '  Error in EDGE.'
    write ( *, '(a,i6)' ) '  IER = ', ier
    stop
  end if
!
!  Test DELNOD by removing nodes 4 to N (in reverse order). 
!
  if ( n <= 3 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Subroutine DELNOD was not tested, because'
    write ( *, '(a)' ) '  the number of nodes N is too small.'

  else

    nn = n

    do

      call delnod ( nn, nn, x, y, z, list, lptr, lend, lnew, nmax, iwk, ier )

      if ( ier /= 0 .and. ier /= 5 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'STRIPACKD_PRB - Fatal error!'
        write ( *, '(a,i6)' ) '  DELNOD returned IER = ', ier
        stop
      end if

      if ( nn <= 3 ) then
        exit
      end if

    end do

  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'STRIPACKD_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  stop
end
