program stripack_prb2
!
!*******************************************************************************
!
!! STRIPACK_PRB2 demonstrates how to use STRIPACK data.
!
!
!  Discussion:
!
!    STRIPACK can compute the Voronoi diagram for data on a sphere.
!
!    This routine has STRIPACK compute the Voronoi diagram, then
!    takes just a few of the "interesting" arrays, and uses them
!    to visit every Voronoi polygon.  Just to prove that is what
!    we are doing, we compute the area of each subtriangle of the
!    polygons, and sum them.  This should be equal to the total area
!    of the sphere, 4 * PI.
!
  implicit none
!
  integer, parameter :: n = 32
!
  double precision d_pi
  integer lend(n)
  integer listc(6*(n-2))
  integer lptr(6*(n-2))
  double precision sphere_area
  double precision x(n)
  double precision xc(2*(n-2))
  double precision y(n)
  double precision yc(2*(n-2))
  double precision z(n)
  double precision zc(2*(n-2))
!
  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'STRIPACK_PRB2'
  write ( *, '(a)' ) '  A sample application of the data produced by'
  write ( *, '(a)' ) '  STRIPACK.  Here, we have STRIPACK compute the'
  write ( *, '(a)' ) '  Voronoi diagram of a set of points on the unit'
  write ( *, '(a)' ) '  sphere, and then we do a simple check by computing'
  write ( *, '(a)' ) '  the polygonal area and summing.'

  call random_number ( harvest = x(1:n) )
  call random_number ( harvest = y(1:n) )
  call random_number ( harvest = z(1:n) )

  call d3vec_normalize ( n, x, y, z )

  call voronoi_get ( n, x, y, z, xc, yc, zc, lend, listc, lptr )

  call voronoi_poly_count ( n, lend, lptr, listc )

  call voronoi_traverse ( n, x, y, z, xc, yc, zc, lend, listc, lptr, &
    sphere_area )

  write ( *, '(a)' ) ' '
  write ( *, '(a,f12.6)' ) &
    '  Sphere area from Voronoi polygons =  ', sphere_area
  write ( *, '(a,f12.6)' ) &
    '  Exact area from spherical geometry = ', 4.0E+00 * d_pi ( )

  stop
end
subroutine voronoi_get ( n, x, y, z, xc, yc, zc, lend, listc, lptr )
!
!*******************************************************************************
!
!! VORONOI_GET calls STRIPACK routines to get Voronoi information.
!
!
!  Modified:
!
!    25 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of points.
!
!    Input, double precision X(N), Y(N), Z(N), the coordinates of points 
!    on the sphere.
!
!    Output, double precision XC(6*(N-2)), YC(6*(N-2)), ZC(6*(N-2)), the 
!    coordinates of the vertices.
!
!    Output, integer LEND(N), points to the "first" vertex in the
!    Voronoi polygon around a particular node.
!
!    Output, integer LPTR(6*(N-2)), given a vertex, returns the next
!    vertex in the Voronoi polygon.  (The vertex numbering is done
!    in such a way that the physical vertex has three distince indices,
!    depending on which polygon we are considering.  Thus, it is always
!    possible to answer the question "which is the next vertex from this
!    one?" because the vertex index also tells us what polygon we are in.)
!
  implicit none
!
  integer n
  integer, parameter :: nrow = 9
!
  double precision ds(n)
  integer i
  integer ierror
  integer iwk(2*n)
  integer lbtri(6,n)
  integer lend(n)
  integer list(6*(n-2))
  integer listc(6*(n-2)) 
  integer lnew
  integer lptr(6*(n-2))
  integer ltri(nrow,2*(n-2))
  integer nb
  double precision norm
  integer nt
  double precision rc(2*(n-2))
  double precision x(n)
  double precision xc(2*(n-2))
  double precision y(n)
  double precision yc(2*(n-2))
  double precision z(n)
  double precision zc(2*(n-2))
!
!  Create the triangulation.
!
  call trmesh ( n, x, y, z, list, lptr, lend, lnew, iwk, iwk(n+1), ds, ierror )

  if ( ierror == -2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'VORONOI_GET - Fatal error!'
    write ( *, '(a)' ) '  Error in TRMESH.'
    write ( *, '(a)' ) '  The first 3 nodes are collinear.'
    stop
  end if

  if ( ierror > 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'VORONOI_GET - Fatal error!'
    write ( *, '(a)' ) '  Error in TRMESH.'
    write ( *, '(a)' ) '  Duplicate nodes encountered.'
    stop
  end if
!
!  Create a triangle list.
!
  call trlist ( n, list, lptr, lend, nrow, nt, ltri, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'VORONOI_GET - Fatal error!'
    write ( *, '(a)' ) '  Error in TRLIST.'
    stop
  end if
!
!  Construct the Voronoi diagram.
!
!  Note that the triangulation data structure is altered if NB > 0.
!
  call crlist ( n, n, x, y, z, list, lend, lptr, lnew, &
    lbtri, listc, nb, xc, yc, zc, rc, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'VORONOI_GET - Fatal error!'
    write ( *, '(a)' ) '  Error in CRLIST.'
    write ( *, '(a,i6)' ) '  IERROR = ', ierror
    stop
  end if

  return
end
subroutine voronoi_traverse ( n, x, y, z, xc, yc, zc, lend, listc, lptr, &
  area_sphere )
!
!*******************************************************************************
!
!! VORONOI_TRAVERSE traverses the polygons in a Voronoi diagram.
!
!
!  Discussion:
!
!    STRIPACK defines a data structure recording the location of
!    the vertices of the Voronoi diagram, and their connectivity.
!    The purpose of this routine is to "visit" each polygon, and,
!    in fact, each subtriangle of each polygon.  Such a procedure
!    would be done when estimating an integral by quadrature, for instance.
!
!  Modified:
!
!    25 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of nodes, and Voronoi polygons.
!
!    Input, double precision X(N), Y(N), Z(N), the coordinates of the nodes.
!
!    Input, double precision XC(6*(N-2)), YC(6*(N-2)), ZC(6*(N-2)), the 
!    coordinates of the vertices.
!
!    Input, integer LEND(N), points to the "first" vertex in the
!    Voronoi polygon around a particular node.
!
!    Input, integer LPTR(6*(N-2)), given a vertex, returns the next
!    vertex in the Voronoi polygon.  (The vertex numbering is done
!    in such a way that the physical vertex has three distince indices,
!    depending on which polygon we are considering.  Thus, it is always
!    possible to answer the question "which is the next vertex from this
!    one?" because the vertex index also tells us what polygon we are in.)
!
!    Output, double precision AREA_SPHERE, the area of the sphere.
!    This is computed simply to demonstrate how the data structure is used.
!    It's assumed the sphere has radius 1, so the area should be
!    4 * PI * R**2 = 12.5664.
!
  implicit none
!
  integer n
!
  double precision area_polygon
  double precision area_sphere
  double precision area_triangle
  double precision areas
  integer lend(n)
  integer listc(6*(n-2))
  integer lptr(6*(n-2))
  integer node
  integer node_last
  integer node_new
  integer node_stop
  double precision r
  double precision v1(3)
  double precision v2(3)
  double precision v3(3)
  integer vertex_last
  integer vertex_new
  double precision x(n)
  double precision xc(6*(n-2))
  double precision y(n)
  double precision yc(6*(n-2))
  double precision z(n)
  double precision zc(6*(n-2))
!
  area_sphere = 0.0E+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'VORONOI_TRAVERSE'
  write ( *, '(a)' ) '  Visit each Voronoi polygon.'
  write ( *, '(a)' ) '  Compute the (spherical) area of each subtriangle'
  write ( *, '(a)' ) '  Add up to get the area of each polygon.'
!
!  To access every polygon, start by accessing a particular node.
!
!  The Voronoi polygon around a node NODE has a pointer LEND(NODE) to the 
!  first (or last) vertex of the Voronoi polygon around NODE.
!
!  To access all the vertices of the polygon in order, start at the
!  special vertex, and then repeatedly use the LPTR array to get the
!  next vertex on the polygon.  Stop when you return to LEND(NODE).
!
!  To subdivide the polygon into triangles, use NODE, VERTEX_LAST,
!  and VERTEX.
!
!  To get the coordinates of these points:
!
!    NODE ==>        X(NODE),         Y(NODE),         Z(NODE).
!
!    VERTEX_LAST ==> XC(VERTEX_LAST), YC(VERTEX_LAST), ZC(VERTEX_LAST)
!    VERTEX      ==> XC(VERTEX     ), YC(VERTEX     ), ZC(VERTEX     ) 
!
  area_sphere = 0.0E+00

  do node = 1, n

    area_polygon = 0.0E+00

    write ( *, '(a)' ) ' '

    node_stop = lend(node)

    node_new = node_stop

    vertex_new = listc(node_new)
!
!  Each iteration of this DO walks along one side of the polygon,
!  considering the subtriangle NODE --> VERTEX_LAST --> VERTEX.
!
    do

      node_last = node_new
      node_new = lptr(node_last)

      vertex_last = vertex_new
      vertex_new = listc(node_new)
!
!  Here is a good place to process information about the polygon side 
!
!   VERTEX_LAST --> VERTEX 
!
!  or about the subtriangle
!
!   NODE --> VERTEX_LAST --> VERTEX.
!
      r = 1.0E+00
      v1(1:3) = (/ x(node),         y(node),         z(node)         /)
      v2(1:3) = (/ xc(vertex_last), yc(vertex_last), zc(vertex_last) /)
      v3(1:3) = (/ xc(vertex_new),  yc(vertex_new),  zc(vertex_new)      /)

      area_triangle = areas ( v1, v2, v3 )

      area_polygon = area_polygon + area_triangle

      write ( *, '(10x,g14.6)' ) area_triangle
!
!  Now if we have reached the vertex where we started, we are done with
!  this polygon.
!
      if ( node_new == node_stop ) then
        exit
      end if

    end do

    write ( *, '(a)' ) ' '
    write ( *, '(i2,g14.6)' ) node, area_polygon

    area_sphere = area_sphere + area_polygon

  end do

  return
end
