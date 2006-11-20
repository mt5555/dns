load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"                      
;
; This code makes an xy plot 
;
begin
; defaults
;
nx = 0
fcor = 0
f_k = 24
epsil_f = .0253303064
;epsil_f = 0.517275
timescale = (epsil_f*(2.*3.141592*f_k)^2)^(.33333)
energyscale = (epsil_f/(2.*3.141592*f_k))^(-2./3.)
;timescale = (epsil_f*(f_k)^2)^(.33333)
;energyscale = (epsil_f/(f_k))^(-2./3.)
print(timescale)
;
;**************************************************
; read in binary data
;**************************************************
  setfileoption("bin","ReadByteOrder","LittleEndian")
;  path  = "../../Projects/KH/Boussinesq/f21norot/hyper4/"
;  path  = "../../Projects/KH/Boussinesq/low/"
  path  = "../../Projects/KH/Boussinesq/B1/"
;  path  = "../../Projects/KH/Boussinesq/n21nf100/"
;  path  = "../../Projects/KH/Boussinesq/r16A1/hyper4/"
;  path  = "../../Projects/KH/Boussinesq/r16A1/"
;  path  = "../../Projects/KH/Boussinesq/n21struct/"
;  path  = "../../Projects/KH/Boussinesq/n21Visc2p5/"
;  path  = "../../Projects/KH/Boussinesq/n21Viscp1/"
;  path  = "../../Projects/KH/Boussinesq/n21Viscp01/"
;  path  = "../../Projects/KH/Boussinesq/n21/"
;  path  = "../../Projects/KH/Boussinesq/n21/hyper4/"
;  path  = "../../Projects/KH/Boussinesq/n21Visc1/"
;  path  = "../../Projects/KH/Boussinesq/n21asp/"
;  path  = "../../Projects/KH/Boussinesq/n21asp/Old/"
;  path  = "../../Projects/KH/Boussinesq/n21asp1/"
;  path  = "../../Projects/KH/Boussinesq/n21aspp1/"
;  path  = "../../Projects/KH/Boussinesq/n21aspp01/"
  pathwrite = path+"ke.txt"
  griddata = cbinread(path+"all.scalars-bous",-1,"double")
  ndata1 = sizeof(griddata)
  ndata2 =ndata1/8. ; 8 bytes per double
  ndata = doubletoint(ndata2)
  print(ndata)
  dims = ndata/10
  print(dims)
  dim = (/dims/)
  time = new(dim,"double")
  ke = new(dim,"double")
  pe = new(dim,"double")
  ke_diss = new(dim,"double")
  pe_diss = new(dim,"double")
  pv = new(dim,"double")
  pv_diss = new(dim,"double")
  potens = new(dim,"double")
  enstr_diss = new(dim,"double")
;  print(griddata)
;
; defaults
; 
nx = 0
fcor = 0
f_k = 0
;
; Loop over all possible data
;
icount = 0
;
;
;
  do while (icount .lt.  dims )
;
;   Assign data
;
    ni = (/griddata(icount*10)/)
    nints = doubletoint(ni)
    time(icount) = (/griddata(icount*10 + 1)/)
    ke(icount) = (/griddata(icount*10 + 2)/)
    pe(icount) = (/griddata(icount*10 + 3)/)
    ke_diss(icount) = (/griddata(icount*10 + 4)/)
    pe_diss(icount) = (/griddata(icount*10 + 5)/)
    pv(icount) = (/griddata(icount*10 + 6)/)
    pv_diss(icount) = (/griddata(icount*10 + 7)/)
    potens(icount) = (/griddata(icount*10 + 8)/)
    enstr_diss(icount) = (/griddata(icount*10 + 9)/)
    icount = icount + 1
;
  end do

  time = time*timescale
  print(time)
  pe = pe*energyscale
  ke = ke*energyscale
  wks = gsn_open_wks("x11","scalars") 
;
; This writes the data out to an asci file. Boy what a pain.
;
;    Define the 2-d matrix
     toprintdims = (/dims,2/)
     toprint = new(toprintdims,"double")
;
;    Write put the kinetic energy
;
     toprint(:,0)=time(:)
     toprint(:,1)=ke(:)
     opt = True
     opt@fout = path+"ke.txt"
     write_matrix(toprint, "2e13.6", opt)
;
;    Write out the potential energy
;
     toprint(:,0)=time(:)
     toprint(:,1)=pe(:)
     opt = True
     opt@fout = path+"pe.txt"
     write_matrix(toprint, "2e13.6", opt)
;
;    Write out the total energy
;
     toprint(:,0)=time(:)
     toprint(:,1)=pe(:) + ke(:)
     opt = True
     opt@fout = path+"tote.txt"
     write_matrix(toprint, "2e13.6", opt)
;
;    Write out the potential vorticity
;
     toprint(:,0)=time(:)
     toprint(:,1)=pv(:)
     opt = True
     opt@fout = path+"pv.txt"
     write_matrix(toprint, "2e13.6", opt)
;
;    Write out the potential enstrophy
;     
     toprint(:,0)=time(:)
     toprint(:,1)=potens(:)
     opt = True
     opt@fout = path+"potens.txt"
     write_matrix(toprint, "2e13.6", opt)
;
;    Write out the potential enstrophy dissipation
;     
     toprint(:,0)=time(:)
     toprint(:,1)=enstr_diss(:)
     opt = True
     opt@fout = path+"potens_diss.txt"
     write_matrix(toprint, "2e15.6", opt)
;
; now plot is ke, pe, tote versus time
;
  y2dims = (/3,dims/)
  print(y2dims)
  y2 = new(y2dims,"double")
  do i = 0,dims-1
   y2(0,i) = ke(i)
   y2(1,i) = pe(i)
   y2(2,i) = ke(i) + pe(i)
  end do
;
  res = True
;  res@xyExplicitLegendLabels = (/"tote","ke","pe"/)
;  res@xyExplicitLabels = (/"A","B","C"/)
  res@tiMainString           = "tote=blue, ke=red, pe=green"
  res@tiXAxisString          = "time"
  res@xyLineColors = (/2,3,4/) ; define line colors
  res@xyLineThicknesses = (/1.,2.,5./) ; define line thickness, 1.0 is deflt
  res@xyDashPattern = (/0.,0.,0./)
  plot = gsn_xy(wks,time,y2,res)

  delete(plot)

end
