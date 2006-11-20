;================================================;
;  conwomap_1.ncl
;================================================;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"   
; ================================================;
begin
;
; Read in the data file
;
  setfileoption("bin","ReadByteOrder","LittleEndian")
  path  = "../../Projects/KH/Boussinesq/B1/"
  udata = cbinread(path+"all.u",-1,"double")
  vdata = cbinread(path+"all.v",-1,"double")
  pvdata = cbinread(path+"all.pv",-1,"double")
;  print(data)
  ndata1 = sizeof(udata)
  ndata2 =ndata1/8. ; 8 bytes per double
  ndata = doubletoint(ndata2)
  print(ndata)
   time = (/udata(0)/)
   print(time)
   nx = doubletoint((/udata(1)/))
   print(nx)
   ny = doubletoint((/udata(2)/))
   print(ny)
   nz = doubletoint((/udata(3)/))
   print(nz)
;
; Declare 3 1-d arrays for the coordinates x,y,z
;
   xcoord = (/udata(4:3+nx)/)
   print(xcoord)
   ycoord = (/udata(4+nx:3+nx+ny)/)
   print(ycoord)
   zcoord = (/udata(4+nx+ny:3+nx+ny+nz)/)
   print(zcoord)
;
; Declare a 2d array for u and v
;
   slicedims = (/nx,ny/)
   uarray = new(slicedims,"double")
   varray = new(slicedims,"double")
   pvarray = new(slicedims,"double")
   print(slicedims)
   icounter = 3+nx+ny+nz
   uarray=0.
   varray=0.
   pvarray = 0.
   do k=1,nz
    do i=1,nx
      do j=1,ny
	icounter = icounter + 1
  	uarray(i-1,j-1) = uarray(i-1,j-1) + (/udata(icounter)/)
  	varray(i-1,j-1) = varray(i-1,j-1) + (/vdata(icounter)/)
  	pvarray(i-1,j-1) = pvarray(i-1,j-1) + (/pvdata(icounter)/)
      end do
    end do
    print(k)
   end do
   do i=1,nx
     do j=1,ny
      uarray(i-1,j-1) = uarray(i-1,j-1)/nz
      varray(i-1,j-1) = varray(i-1,j-1)/nz
      pvarray(i-1,j-1) = pvarray(i-1,j-1)/nz
     end do
   end do
;
; Min Max of the data 
;
  data_max = max(pvarray)
  data_min = min(pvarray)
  print(data_max)
  print(data_min)
;=================================================;
; PLOT 2
;=================================================;
  wks  = gsn_open_wks("X11","vector")     ; open a ps file
;  wks  = gsn_open_wks("eps","vector")     ; open a ps file
  setvalues NhlGetWorkspaceObjectId()
  "wsMaximumSize": 33554432
  end setvalues  
  gsn_define_colormap(wks,"gui_default")  ; choose a colormap

  res                     = True          ; plot mods desired

  res@cnFillOn            = True          ; turn on color
  res@lbTitleFontHeightF  = .007              ; controls the font size 
  res@lbOrientation =   "vertical"          ; vertical label bar
  res@gsnSpreadColors     = True          ; use full range of colormap
 
;  res@vpWidthF            = 0.5           ; change aspect ratio of plot
;  res@vpHeightF           = 0.3

  plot=gsn_csm_vector_scalar_map_ce(wks,uarray,varray,\
                                    pvarray,res)

; note that the default plot draw here still uses the long_name and
; units of the variable as labels.

; In this file, the coordinate variables of the data are listed as 
; lat and long, but they are really just index points, which makes this
; data suitable for plotting without a map.
  delete(wks)
  end
