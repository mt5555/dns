load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"                      
;
; This code plots all the spectra in a .spec like file 
;
begin
;
; Read in the kinetic energy spectrum
;
  setfileoption("bin","ReadByteOrder","LittleEndian")
;  path  = "../../Projects/KH/Boussinesq/n21/"
;  path  = "../../Projects/KH/Boussinesq/n21/hyper4/"
;  path  = "../../Projects/KH/Boussinesq/r16A1/hyper4/"
;  path  = "../../Projects/KH/Boussinesq/B1/"
; path  = "../../Projects/KH/Boussinesq/low/"
  path  = "../../Projects/Boussinesq/kd/kd40/Data/"
  data = cbinread(path+"all.spec",-1,"double")
;  print(data)
;
; Now read in the potential energy spectrum
;
  datap = cbinread(path+"all.pspec",-1,"double")
;;  print(datap)
;
; Now examine the kinetic energy spectrum data
;
  ndata1 = sizeof(data)
  ndata2 =ndata1/8. ; 8 bytes per double
  ndata = doubletoint(ndata2)
  print(ndata)
  time =  (/data(0)/)
  print(time)
  nr =  (/data(1)/)
  print(nr)
  nri = doubletoint(nr)
  nx = (/data(nri+2)/)
  nxi = doubletoint(nx)
  print(nxi)
  ny = (/data(nri+3*nxi+3)/)
  nyi = doubletoint(ny)
  print(nyi)
  nz = (/data(nri+3*nxi+3*nyi+4)/)
  nzi = doubletoint(nz)
  print(nzi)
  ndatapertime = nri+5+3*nxi+3*nyi+3*nzi
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
   specparraym1 = new(rdimsm1,"double")
   totenarraym1 = new(rdimsm1,"double")
;
;  Assign the kinetic energy spectrums to other arrays
;

  nskp = 0
  iterindex = 0
  do while (nskp .lt. ndata) 
     
     iterindex = iterindex + 1
 ;;    print(iterindex)
     timedata(iterindex-1) =  (/data(nskp + 0)/)
     nr =  (/data(nskp+1)/)
     spec_r = (/data(nskp +2:nskp + nri + 1 )/)
     nx = (/data(nskp +nri+2)/)
     nxi = doubletoint(nx)
     specux=(/data(nskp+3+nri:nskp+3+nri+nxi)/)
     specuy=(/data(nskp+3+nri+nxi:nskp+2+nri+2*nxi)/)
     specuz=(/data(nskp+3+nri+2*nxi:nskp+2+nri+3*nxi)/)
     ny = (/data(nskp+nri+3*nxi)/)
     nyi = doubletoint(ny)
     specvx=(/data(nskp+4+nri+3*nxi:nskp+3+nri+3*nxi+nyi)/)
     specvy=(/data(nskp+4+nri+3*nxi+nyi:nskp+3+nri+3*nxi+2*nyi)/)
     specvz=(/data(nskp+4+nri+3*nxi+2*nyi:nskp+3+nri+3*nxi+3*nyi)/)
     nz = (/data(nskp+nri+3*nxi+3*nyi)/)
     nzi = doubletoint(nz)
     specvx=(/data(nskp+5+nri+3*nxi+3*nyi:nskp+4+nri+3*nxi+3*nyi+nzi)/)
     specvy=(/data(nskp+5+nri+3*nxi+3*nyi+nzi:nskp+4+nri+3*nxi+3*nyi+2*nzi)/)
     specvz=(/data(nskp+5+nri+3*nxi+3*nyi+2*nzi:nskp+4+nri+3*nxi+3*nyi+3*nzi)/)
     nskp = nskp + ndatapertime
; 
;    Since I suck at assigning values for multi-d arrays I brute force
;    it here. I'm assigning it to a 2-d array here for plotting
;     multiple curves, but I may not have to do it this way if I 
;    use overlay correctly.
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
    opt@fout = path+timedata(iterindex-1)+"kespec.txt"
    write_matrix(toprint, "2e13.6",opt)
	

  end do
;
; Redistribute the potential energy
;


  ndatap1 = sizeof(datap)
  ndatap2 =ndatap1/8. ; 8 bytes per double
  ndatap = doubletoint(ndatap2)
  print(ndatap)
  nscalars = (/datap(0)/)
  time =  (/datap(1)/)
  nrp =  (/datap(2)/)
  nrpi = doubletoint(nrp)
;;  print(nrpi)
  nxp = (/datap(nrpi+3)/)
  nxpi = doubletoint(nxp)
;;  print(nxpi)
  nyp = (/datap(nrpi+nxpi+4)/)
  nypi = doubletoint(nyp)
;;  print(nypi)
  nzp = (/datap(nrpi+nxpi+nypi+5)/)
  nzpi = doubletoint(nzp)
;;  print(nzpi)
  ndatapertimep = nrpi+6+nxpi+nypi+nzpi
  print(ndatapertimep)
  ntotaldatap = ndatap/ndatapertimep ; This is the number of time
                                     ; snapshots that were read in
  ntp = doubletoint(ntotaldatap)
  print(ntotaldatap)

   rdims = (/ntp,nrpi/)
   xdims = (/ntp,nxpi/)
   ydims = (/ntp,nypi/)
   zdims = (/ntp,nzpi/)
   specprarray = new(rdims,"double")
   totenarray = new(rdims,"double")
;
; Assign all the theta data to arrays
;
  nskp = 0
  iterindex = 0
  do while (nskp .lt. ndatap) 
     
     iterindex = iterindex + 1
     nscalars = (/datap(nskp + 0)/)
     timedata(iterindex-1) =  (/datap(nskp + 1)/)
     nrp =  (/datap(nskp+2)/)
     nrpi = doubletoint(nrp)
     specp_r = (/datap(nskp + 3:nskp + nrpi + 2 )/)
     nxp = (/datap(nskp+nrpi+3)/)
     nxpi = doubletoint(nxp) 
     specpx=(/datap(nskp+4+nrpi:nskp+3+nrpi+nxpi)/)
     nyp = (/datap(nskp+4+nrpi+nxpi)/)
     nypi = doubletoint(nyp)
     specpy=(/datap(nskp+5+nrpi+nxpi:nskp+4+nrpi+nxpi+nypi)/)
     nzp = (/datap(nskp+5+nrpi+nxpi+nypi)/)
     nzpi = doubletoint(nzp)
     specpz=(/datap(nskp+6+nrpi+nxpi+nypi:nskp+5+nrpi+nxpi+nypi+nzpi)/)
     nskp = nskp + ndatapertimep
; 
;    Since I suck at assigning values for multi-d arrays I brute force
;    it here. I'm assigning it to a 2-d array here for plotting
;     multiple curves, but I may not have to do it this way if I 
;    use overlay correctly.
;
     do i = 1,nri
      specprarray(iterindex-1,i-1) = (/specp_r(i-1)/)
      totenarray(iterindex-1,i-1) = specprarray(iterindex-1,i-1) + specrarray(iterindex-1,i-1)
     end do
     if (iterindex.gt.1) then
       do i = 1,nri
         specparraym1(iterindex-2,i-1) = (/specp_r(i-1)/)
         totenarraym1(iterindex-2,i-1) = specarraym1(iterindex-2,i-1) + (/specp_r(i-1)/)
       end do
     end if
;
;   Write the pe spec to asci files
;
    toprint(:,0)=krnumbers(:)
    toprint(:,1)=specp_r(:)
    opt=True
    opt@fout = path+timedata(iterindex-1)+"pespec.txt"
    write_matrix(toprint, "2e13.6",opt)
;
;   Write the totenergy spec to asci files
;
    toprint(:,0)=krnumbers(:)
    toprint(:,1)= specp_r(:)+spec_r(:)
    opt=True
    opt@fout = path+timedata(iterindex-1)+"totespec.txt"
    write_matrix(toprint, "2e13.6",opt)

  end do
;;;
;;; Plot the total energy spectrum
;;;
;
; Open a workstation
;
  wks = gsn_open_wks("x11","totspectrum") 
;
; Set some resources
;
  res = True
  res@tiMainString = "Total energy spectra"  
  res@tiXAxisString = "k" 
  res@tiYAxisString = "Etotal" 
  res@xyYStyle = "Log"
  res@xyXStyle = "Log"
;;  res@xyComputeYMin = "True"
  res@trYMinF = .00000001
  res@trYMaxF = 2.0
  res@trXMinF = 1
;
; Plot the most recent
;
  plot = gsn_csm_xy(wks,krnumbers,totenarray(ntp-1,:),res)
;
; Plot them all
;
  plot = gsn_csm_xy(wks,krnumbers,totenarraym1,res)
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
  res@tiMainString = "Total Energy Spectrum"  
  res@gsnDraw      = False                        ; don't draw yet
  res@gsnFrame     = False                        ; don't advance frame yet
  res@xyLineColors = (/"blue"/) ; define line colors
  firstplot = gsn_csm_xy(wks,krnumbers,totenarray(1,:),res)
  res@xyLineColors = (/"green"/) ; define line colors
  res@pmLegendParallelPosF = .8 ; move units right
  res@pmLegendOrthogonalPosF = -.6 ; move units down
  res@xyExplicitLegendLabels = (/"first"/)
  secondplot = gsn_csm_xy(wks,krnumbers,totenarray(ntp-1,:),res)
  overlay(firstplot,secondplot)
  draw(firstplot)
  frame(wks)
  delete(plot)
  delete(firstplot)
  delete(secondplot)
;;;
;;; Plot the kinetic energy spectrum
;;;
;
; Set some resources
;
  reske = True
  reske@tiMainString = "Kinetic energy spectra"  
  reske@tiXAxisString = "k" 
  reske@tiYAxisString = "KE" 
  reske@xyYStyle = "Log"
  reske@xyXStyle = "Log"
;;  res@xyComputeYMin = "True"
  reske@trYMinF = .00000001
  reske@trYMaxF = 2.0
  reske@trXMinF = 1
;
; Plot the most recent
;
  plotke = gsn_csm_xy(wks,krnumbers,specarraym1(ntp-2,:),reske)
;
; Plot them all
;
  plotke = gsn_csm_xy(wks,krnumbers,specarraym1,reske)
;
; Just overlay the first and the last spectra
;
  reske@pmLegendDisplayMode = "Always"              ;turn on the legend
  reske@pmLegendSide = "BottomRight"
  reske@pmLegendWidthF=.1
  reske@pmLegendHeightF=.4
  reske@lgLabelFontHeightF = .02
  reske@pmLegendSide="Top"
  reske@pmLegendParallelPosF = .8 ; move units right
  reske@pmLegendOrthogonalPosF = -.5 ; move units down
  reske@lgPerimOn = False ; turn off the box around the legend
  reske@xyExplicitLegendLabels = (/"recent"/)
  reske@trYMaxF = 1.
  reske@tiMainString = "Kinetic Energy Spectrum"  
  reske@gsnDraw      = False                        ; don't draw yet
  reske@gsnFrame     = False                        ; don't advance frame yet
  reske@xyLineColors = (/"blue"/) ; define line colors
  firstplotke = gsn_csm_xy(wks,krnumbers,specarraym1(1,:),reske)
  reske@xyLineColors = (/"green"/) ; define line colors
  reske@pmLegendParallelPosF = .8 ; move units right
  reske@pmLegendOrthogonalPosF = -.6 ; move units down
  reske@xyExplicitLegendLabels = (/"first"/)
  secondplotke = gsn_csm_xy(wks,krnumbers,specarraym1(ntp-2,:),reske)
  overlay(firstplotke,secondplotke)
  draw(firstplotke)
  frame(wks)
  delete(plotke)
  delete(firstplotke)
  delete(secondplotke)
;;;
;;; Plot the potential energy spectrum
;;;
;
; Set some resources
;
  respe = True
  respe@tiMainString = "Potential energy spectra"  
  respe@tiXAxisString = "k" 
  respe@tiYAxisString = "PE" 
  respe@xyYStyle = "Log"
  respe@xyXStyle = "Log"
;;  res@xyComputeYMin = "True"
  respe@trYMinF = .00000001
  respe@trYMaxF = 2.0
  respe@trXMinF = 1
;
; Plot the most recent
;
  plotpe = gsn_csm_xy(wks,krnumbers,specparraym1(ntp-2,:),respe)
;
; Plot them all
;
  plotpe = gsn_csm_xy(wks,krnumbers,specparraym1,respe)
;
; Just overlay the first and the last spectra
;
  respe@pmLegendDisplayMode = "Always"              ;turn on the legend
  respe@pmLegendSide = "BottomRight"
  respe@pmLegendWidthF=.1
  respe@pmLegendHeightF=.4
  respe@lgLabelFontHeightF = .02
  respe@pmLegendSide="Top"
  respe@pmLegendParallelPosF = .8 ; move units right
  respe@pmLegendOrthogonalPosF = -.5 ; move units down
  respe@lgPerimOn = False ; turn off the box around the legend
  respe@xyExplicitLegendLabels = (/"recent"/)
  respe@trYMaxF = 1.
  respe@tiMainString = "Potential Energy Spectrum"  
  respe@gsnDraw      = False                        ; don't draw yet
  respe@gsnFrame     = False                        ; don't advance frame yet
  respe@xyLineColors = (/"blue"/) ; define line colors
  firstplotpe = gsn_csm_xy(wks,krnumbers,specarraym1(1,:),respe)
  respe@xyLineColors = (/"green"/) ; define line colors
  respe@pmLegendParallelPosF = .8 ; move units right
  respe@pmLegendOrthogonalPosF = -.6 ; move units down
  respe@xyExplicitLegendLabels = (/"first"/)
  secondplotpe = gsn_csm_xy(wks,krnumbers,specparraym1(ntp-2,:),respe)
  overlay(firstplotpe,secondplotpe)
  draw(firstplotpe)
  frame(wks)
  delete(plotpe)
  delete(firstplotpe)
  delete(secondplotpe)

end
