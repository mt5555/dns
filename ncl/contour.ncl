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
  setfileoption("bin","ReadByteOrder","LittleEndian")
  path  = "../../Projects/KH/Boussinesq/B1/"
  data = cbinread(path+"all.pv",-1,"double")
;  print(data)
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
   print(xcoord)
   ycoord = (/data(4+nx:3+nx+ny)/)
   print(ycoord)
   zcoord = (/data(4+nx+ny:3+nx+ny+nz)/)
   print(zcoord)
;
; Declare a 2d array to test
;
   
   slicedims = (/nx,ny/)
   contarray = new(slicedims,"double")
   print(slicedims)
   icounter = 3+nx+ny+nz
   contarray=0
   navg = nz
   do k=1,navg
    do i=1,nx
      do j=1,ny
	icounter = icounter + 1
  	contarray(i-1,j-1) = contarray(i-1,j-1) + (/data(icounter)/)
      end do
    end do
    print(k)
   end do
   do i=1,nx
     do j=1,ny
      contarray(i-1,j-1) = contarray(i-1,j-1)/navg
;      if (fabs((/contarray(i-1,j-1)/)).gt.1.) then
;         print(i)
;	 print(j)
;         print((/contarray(i-1,j-1)/))
;      end if
     end do
   end do
;
; Min Max of the data 
;
  data_max = max(contarray)
  data_min = min(contarray)
  print(data_max)
  print(data_min)
;=================================================;
; PLOT 2
;=================================================;
;  wks  = gsn_open_wks("X11","cont")     ; open a ps file
  wks  = gsn_open_wks("eps","cont")     ; open a ps file
  setvalues NhlGetWorkspaceObjectId()
  "wsMaximumSize": 33554432
  end setvalues  
  gsn_define_colormap(wks,"gui_default")  ; choose a colormap

  res                     = True          ; plot mods desired

;  The following is for setting manual levels
;  res@cnLevelSelectionMode = "ManualLevels"
;  res@cnLevelSpacingF     = .1
;  res@cnMinLevelValF      = -.4
;  res@cnMaxLevelValF      = .4
  res@cnFillOn            = True          ; turn on color
  res@lbTitleFontHeightF  = .007              ; controls the font size 
  res@lbOrientation =   "vertical"          ; vertical label bar
  res@gsnSpreadColors     = True          ; use full range of colormap
  res@cnLinesOn          = False 

;  res@vpWidthF            = 0.5           ; change aspect ratio of plot
;  res@vpHeightF           = 0.3

  plot = gsn_csm_contour(wks,contarray,res)      ; contour the variable

; note that the default plot draw here still uses the long_name and
; units of the variable as labels.

; In this file, the coordinate variables of the data are listed as 
; lat and long, but they are really just index points, which makes this
; data suitable for plotting without a map.
  delete(wks)
  end
