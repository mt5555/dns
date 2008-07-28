load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"                      
;
; This code plots all the spectra in a .pv2spec like file 
;
begin
  setfileoption("bin","ReadByteOrder","LittleEndian")
  path  = "../../Projects/KH/Boussinesq/np08rp08/"
  path  = "../../Projects/Boussinesq/kd/kd40/Data/"
  data = cbinread(path+"all.pv2spec",-1,"double")
;  print(data)
  ndata1 = sizeof(data)
  ndata2 =ndata1/8. ; 8 bytes per double
  ndata = doubletoint(ndata2)
;  print(ndata)
  time =  (/data(0)/)
;  print(time)
  nr =  (/data(1)/)
  nri = doubletoint(nr)
;  print(nri)
  nx = (/data(2)/)
  nxi = doubletoint(nx)
;  print(nxi)
  ny = (/data(3)/)
  nyi = doubletoint(ny)
;  print(nyi)
  nz = (/data(4)/)
  nzi = doubletoint(nz)
;  print(nzi)
  iwave2d = (/data(5)/)
  iwave2di = doubletoint(iwave2d)
;  print(iwave2di)
  ndatapertime = nri+6+nxi+nyi+nzi+iwave2di*nzi
  print(ndatapertime)
  ntotaldata = ndata/ndatapertime ; This is the number of time
                                  ;snapshots that were read in
  nt = doubletoint(ntotaldata)
  print(ntotaldata)
;
;  Allocate some memory
;
   timedata = new(nt,"double")
   rdims = (/nt,nri/)
   xdims = (/nt,nxi/)
   ydims = (/nt,nyi/)
   zdims = (/nt,nzi/)
   specrarray = new(rdims,"double")
;
; Define the k grid for the spherical shells
;
  dkr = 1 
  krmin = 1
  krnumbers = krmin + ispan(0,nri-1,1)*dkr
;
; Do the same for x, y, and z
;
  dkx = 1
  kxmin = 1
  kxnumbers = kxmin + ispan(0,nxi-1,1)*dkx
  dky = 1
  kymin = 1
  kynumbers = kymin + ispan(0,nyi-1,1)*dky
  dkz = 1
  kzmin = 1
  kznumbers = kzmin + ispan(0,nzi-1,1)*dkz
;
; Make an array for printing ascii files
;
   toprintdims = (/nri,2/)
   toprint = new(toprintdims,"double")
;
;  Also make an array that will exclude the initial spectra
;  since its zero and causes trouble with the plotting
   rdimsm1 = (/nt-1,nri/)
   specarraym1 = new(rdimsm1,"double")
   spec_r = (/data(6:nri+5)/)
   specux = (/data(6+nri:5+nri+nxi)/)
   specuy = (/data(6+nri+nxi:5+nri+nxi+nyi)/)
   specuz = (/data(6+nri+nxi+nyi:5+nri+nxi+nyi+nzi)/)
   do i = 1,nzi
      spec_2d = (/data(6+i-1+nri+nxi+nyi+nzi+(i-1)*iwave2di:5+i-1+nri+nxi+nyi+nzi+i*iwave2di)/)
   end do


   nskp = 0
   iterindex = 0
   do while (nskp .lt. ndata) 
     
     iterindex = iterindex + 1
;     print(iterindex)
     timedata(iterindex-1) =  (/data(nskp + 0)/)
     print(timedata(iterindex-1))
     nr =  (/data(nskp+1)/)
     nri = doubletoint(nr)
;     print(nri)
     nx = (/data(nskp +2)/)
     nxi = doubletoint(nx)
     print(nxi)
     ny = (/data(nskp+3)/)
     nyi = doubletoint(ny)
     print(nyi)
     nz = (/data(nskp+4)/)
     nzi = doubletoint(nz)
     print(nzi)
     iwave2d = (/data(nskp+5)/)
     iwave2di = doubletoint(iwave2di)
     print(iwave2di)
     spec_r = (/data(nskp+6:nskp+nri+5)/)
     print(spec_r(0))
     specux = (/data(nskp+6+nri:nskp+5+nri+nxi)/)
     specuy = (/data(nskp+6+nri+nxi:nskp+5+nri+nxi+nyi)/)
     specuz = (/data(nskp+6+nri+nxi+nyi:nskp+5+nri+nxi+nyi+nzi)/)
     do i = 1,nzi
        spec_2d = (/data(nskp+6+nri+nxi+nyi+nzi+(i-1)*iwave2di:nskp+5+nri+nxi+nyi+nzi+i*iwave2di )/)
     end do
;     print(spec_2d(iwave2di-1))

     nskp = nskp + ndatapertime

; 
;    Since I suck at assigning values for multi-d arrays I brute force
;    it here.
;
     do i = 1,nri
      specrarray(iterindex-1,i-1) = (/spec_r(i-1)/)
     end do
     if (iterindex.gt.1) then
       do i = 1,nri
         specarraym1(iterindex-2,i-1) = (/spec_r(i-1)/)
       end do
     end if

;
;   Write the data out to asci files
;
    toprint(:,0)=krnumbers(:)
    toprint(:,1)=spec_r(:)
    opt=True
    opt@fout = path+timedata(iterindex-1)+"pv2spec.txt"
    write_matrix(toprint, "2e13.6",opt)

  end do
;
; Open a workstation
;

  wks = gsn_open_wks("x11","spectrum") 
;
; Set some resources
;

  res = True
  res@tiMainString = "Potential Enstrophy spectra"  
  res@tiXAxisString = "k" 
  res@tiYAxisString = "Q2" 
  res@xyYStyle = "Log"
  res@xyXStyle = "Log"
  res@trYMinF = .0000001
  res@trXMinF = 1
;
; Plot the most recent
;
  res@tiMainString = "Potential Enstrophy spectra - recent time"  
  plot = gsn_csm_xy(wks,krnumbers,spec_r,res)
;
; Plot them all
;
  res@tiMainString = "Potential Enstrophy spectra - all times"  
  plot = gsn_csm_xy(wks,krnumbers,specarraym1,res)
;
; Just overlay the first and the last spectra
;
  res@pmLegendDisplayMode = "Always"              ;turn on the legend
  res@pmLegendSide = "BottomRight"
  res@pmLegendWidthF=.1
  res@pmLegendHeightF=.4
  res@lgLabelFontHeightF = .02
  res@pmLegendSide="Top"
  res@pmLegendParallelPosF = .8 ; move units right
  res@pmLegendOrthogonalPosF = -.5 ; move units down
  res@lgPerimOn = False ; turn off the box around the legend
   res@xyExplicitLegendLabels = (/"recent"/)
  res@trYMaxF = 1.
  res@tiMainString = "Potential Enstrophy spectra - first and last"  
  res@gsnDraw      = False                        ; don't draw yet
  res@gsnFrame     = False                        ; don't advance frame yet
  res@xyLineColors = (/"blue"/) ; define line colors
  firstplot = gsn_csm_xy(wks,krnumbers,spec_r,res)
  res@xyLineColors = (/"green"/) ; define line colors
  res@pmLegendParallelPosF = .8 ; move units right
  res@pmLegendOrthogonalPosF = -.6 ; move units down
  res@xyExplicitLegendLabels = (/"first"/)
  secondplot = gsn_csm_xy(wks,krnumbers,specarraym1(0,:),res)
  overlay(firstplot,secondplot)
  draw(firstplot)
  frame(wks)

  delete(plot)


end
