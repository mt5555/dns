;================================================;
;  contour.ncl
;================================================;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"   
; ================================================;
begin
;
; Read in the data file
;
; ivert=0 means horizontal slices, ivert=1 means vertical slices
  ivert = 0
  nhslice = 128
  nvslice = 128
; ipotens = 0 means don't compute the potential enstrophy from the pv
  ipotens = 0
;
; Read in the data
;
  setfileoption("bin","ReadByteOrder","LittleEndian")
  path  = "../../Projects/KH/Boussinesq/128/n100f100/"
  data = cbinread(path+"n100f1000004.0000.u",-1,"double")
  ndata1 = sizeof(data)
  ndata2 =ndata1/8. ; 8 bytes per double
  ndata = doubletoint(ndata2)
  print(ndata)
   time = (/data(0)/)
   print(time)
   nx = doubletoint((/data(1)/))
   print(nx)
   ny = doubletoint((/data(2)/))
   print(ny)
   nz = doubletoint((/data(3)/))
   print(nz)
;
; Declare 3 1-d arrays for the coordinates x,y,z
;
   xcoord = (/data(4:3+nx)/)
;   print(xcoord)
   ycoord = (/data(4+nx:3+nx+ny)/)
;   print(ycoord)
   zcoord = (/data(4+nx+ny:3+nx+ny+nz)/)
;   print(zcoord)
;
; Declare a 3d array for storing the data
;
   bigarraydims = (/nx,ny,nz/)
   bigarray = new(bigarraydims,"double")
   average_value = new(1,"double")
   max_value = new(1,"double")
   abs_value = new(1,"double")
;
;  Store the 3d data in a way we can do something with.
;
   bigarray=0
;
   icounter = 3+nx+ny+nz
   do k=1,nz
    do i=1,nx
      do j=1,ny
	icounter = icounter + 1
           bigarray(i-1,j-1,k-1) = (/data(icounter)/)
       end do
    end do
   end do
;
;    Find averages
;
     max_value = max(fabs(bigarray))
     average_value = avg(bigarray)
     abs_value = avg(fabs(bigarray))
     print(average_value)
     print(abs_value)
     print(max_value)

  end
