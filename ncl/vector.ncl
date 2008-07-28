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
  nslice = 64
  setfileoption("bin","ReadByteOrder","LittleEndian")
  path  = "../../Projects/Boussinesq/kd/kd40/Data/"
  udata = cbinread(path+"kd400032.0000.u",-1,"double")
  vdata = cbinread(path+"kd400032.0000.v",-1,"double")
  pvdata = cbinread(path+"kd400032.0000.pv",-1,"double")
;  print(data)
  ndata1 = sizeof(udata)
  ndata2 =ndata1/8. ; 8 bytes per double
  ndata = doubletoint(ndata2)
;  print(ndata)
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
; Declare a 3d array for u and v and potens
;
   slicedims = (/nx,ny,nz/)
   uarray = new(slicedims,"double")
   varray = new(slicedims,"double")
   pvarray = new(slicedims,"double")
   print(slicedims)
   icounter = 3+nx+ny+nz
   uarray=0.
   varray=0.
   pvarray = 0.
   navg = 2
   do k=1,nz
     do i=1,nx   
       do j=1,ny
	 icounter = icounter + 1
   	 uarray(i-1,j-1,k-1) = (/udata(icounter)/)
  	 varray(i-1,j-1,k-1) = (/vdata(icounter)/)
  	 pvarray(i-1,j-1,k-1) = ((/pvdata(icounter)/)^2)/2
       end do
     end do
;     print(k)
   end do
;
;  Reduce the number of grid points so as to see the vectors
;
  nx_small = (nx-1)/4
  ny_small = (ny-1)/4
  nz_small = (nz-1)/4
  nslice_small = nz_small/2
;
; Declare a 2d array for u and v and potens
;
   slicedims2 = (/nx_small,ny_small/)
   print(slicedims2)
   u2d = new(slicedims2,"double")
   v2d = new(slicedims2,"double")
   potens2d = new(slicedims2,"double")

   do i=1,nx
     do j = 1,ny
       u2d = uarray(i-1,j-1,nslice-1)
       v2d = varray(i-1,j-1,nslice-1)
       potens2d = pvarray(i-1,j-1,nslice-1)
     end do
   end do
;
; Min Max of the 2d data 
;
  u_max = max(u2d)
  u_min = min(u2d)
  print(u_max)
  print(u_min)
  v_max = max(v2d)
  v_min = min(v2d)
  print(v_max)
  print(v_min)
  potens_max = max(pvarray)
  potens_min = min(pvarray)
  print(potens_max)
  print(potens_min)
;
;=================================================;
; PLOT 2
;=================================================;
;  wks  = gsn_open_wks("X11","vector")     ; open a ps file
  wks  = gsn_open_wks("eps","vector")     ; open a ps file
  setvalues NhlGetWorkspaceObjectId()
  "wsMaximumSize": 33554432
  end setvalues  
  gsn_define_colormap(wks,"gui_default")  ; choose a colormap

  res                     = True          ; plot mods desired

  res@lbTitleFontHeightF  = .007              ; controls the font size 
  res@lbOrientation =   "vertical"          ; vertical label bar
  res@gsnSpreadColors     = True          ; use full range of colormap
;  res@vcMinFracLengthF     = .30
;  res@vcRefMagnitudeF     = .10
;  res@vcRefLengthF        = .30
;  res@vcLevelSelectionMode = "ManualLevels"
  res@vcLevelSelectionMode = "AutomaticLevels"
;  res@vcLevelSpacingF      = .20
;  res@vcMinLevelValF       = -.5
;  res@vcMaxLevelValF       = .5
;   res@vcFillArrowsOn           = True   
;   res@vcMonoFillArrowFillColor = False ; in different colors
;   res@vcFillArrowEdgeColor     = 1     ; Draw the edges in black.
;   res@vcFillArrowWidthF        = 0.055 ; Make vectors thinner.
  res@vpWidthF            = 2.0           ; change aspect ratio of plot
  res@vpHeightF           = 2.0
  plot1=gsn_vector(wks,u2d,v2d,res)
;  plot=gsn_csm_vector_scalar(wks,uarray,varray,\
;                                    pvarray,res)

  delete(wks)
  end
