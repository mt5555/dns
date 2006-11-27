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
  path  = "../../Projects/KH/Boussinesq/B1/"
  data = cbinread(path+"B1.pv",-1,"double")
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
;
; Declare a 2d array for contouring and analysis
;
   slicedims = (/nx,ny/)
   contarray = new(slicedims,"double")
   averagearray = new(slicedims,"double")
   diffarray = new(slicedims,"double")
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
   print(ivert)
;
  if (ivert .eq. 0) then
;
;    The following is for doing vertical averages
;
     contarray=0
     averagearray=0
     navg = nz
     do k=1,navg
      do i=1,nx
        do j=1,ny
        if (ipotens.eq.0) then
    	  averagearray(i-1,j-1) = averagearray(i-1,j-1) + bigarray(i-1,j-1,k-1)
        else
    	  averagearray(i-1,j-1) = averagearray(i-1,j-1) + bigarray(i-1,j-1,k-1)*bigarray(i-1,j-1,k-1)/2.
        end if          
        end do
      end do
     end do
;
;    Normalize
;
     do i=1,nx
       do j=1,ny
        averagearray(i-1,j-1) = averagearray(i-1,j-1)/navg
       end do
     end do
     print(min(averagearray))
     print(max(averagearray))
;
;    Pick out a slice
;
     contarray = 0.
      do i=1,nx
        do j=1,ny
        if (ipotens.eq.0) then
    	  contarray(i-1,j-1) = bigarray(i-1,j-1,nhslice-1)
        else
    	  contarray(i-1,j-1) = bigarray(i-1,j-1,nhslice-1)*bigarray(i-1,j-1,nhslice-1)/2.
        end if
        end do
      end do
     print(min(contarray))
     print(max(contarray))
;
;    See how much its different from the average
;
      do i=1,nx
        do j=1,ny
    	  diffarray(i-1,j-1) = contarray(i-1,j-1) - averagearray(i-1,j-1)
        end do
      end do
     print(min(diffarray))
     print(max(diffarray))
;
;
;
   else
;
;  Look at vertical slices
;
     averagearray=0
     contarray=0
;
;    Compute the horizontal average
;
     do i=1,nx
       do k=1,nz
         do j=1,ny
            if(ipotens.eq.0) then
              averagearray(j-1,k-1) = averagearray(j-1,k-1) + bigarray(i-1,j-1,k-1)
            else
              averagearray(j-1,k-1) = averagearray(j-1,k-1) + bigarray(i-1,j-1,k-1)*bigarray(i-1,j-1,k-1)/2.
            end if
          end do
        end do
      end do
;
;     Normalize
;
      do k=1,nz
        do j=1,ny
         averagearray(j-1,k-1) = averagearray(j-1,k-1)/nx
        end do
      end do
;
;     Pick out a slice to look at
;
      do k=1,nz
        do j=1,ny
        if (ipotens.eq.0) then
         contarray(j-1,k-1) = bigarray(nvslice-1,j-1,k-1)
        else
         contarray(j-1,k-1) = bigarray(nvslice-1,j-1,k-1)*bigarray(nvslice-1,j-1,k-1)/2.
        end if
        end do
      end do
; 
;     Find out how far your snapshot is from the average
;
      diffarray=0
      do k=1,nz
      do j=1,ny
          diffarray(j-1,k-1) = contarray(j-1,k-1)-averagearray(j-1,k-1)
      end do
      end do
;
     end if
;
; Min Max of the data 
;
  print(max(contarray))
  print(min(contarray))
  print(max(averagearray))
  print(min(averagearray))
  print(max(diffarray))
  print(min(diffarray))
;=================================================;
; PLOT 1
;=================================================;
;  wks  = gsn_open_wks("X11","cont")     ; open a ps file
  wks  = gsn_open_wks("eps","cont")     ; open a ps file for the snapshot
  wks2  = gsn_open_wks("eps","contavg")     ; open a ps file for the average
  wks3  = gsn_open_wks("eps","contdiff")     ; open a ps file for the difference
  setvalues NhlGetWorkspaceObjectId()
  "wsMaximumSize": 500000000
  end setvalues  
  gsn_define_colormap(wks,"gui_default")  ; choose a colormap
  gsn_define_colormap(wks2,"gui_default")  ; choose a colormap
  gsn_define_colormap(wks3,"gui_default")  ; choose a colormap

  res                     = True          ; plot mods desired

;  res@cnLevelSelectMode   = "ExplicitLevels"
;  res@cnLevels            = (/-3.,-2.,-1.,1.,2.,3./)
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

  plot = gsn_csm_contour(wks,contarray,res)      ; contour the snapshot
  delete(wks)
  plot2 = gsn_csm_contour(wks2,averagearray,res)      ; contour the average
  delete(wks2)
  plot3 = gsn_csm_contour(wks3,diffarray,res)      ; contour the difference
  delete(wks3)
;
; note that the default plot draw here still uses the long_name and
; units of the variable as labels.

; In this file, the coordinate variables of the data are listed as 
; lat and long, but they are really just index points, which makes this
; data suitable for plotting without a map.
  end
