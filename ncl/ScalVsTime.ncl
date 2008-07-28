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
f_k = 3
;epsil_f = .0253303064
epsil_f = 0.517275
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
;  path  = "../../Projects/KH/Boussinesq/Bp01/"
;  path  = "../../Projects/Boussinesq/512/rop005fr1/Data/"
  path  = "../../Projects/Boussinesq/640/n1f1000/"
  path  = "../../Projects/Boussinesq/128/rop05fr1kf2/Scalars/"
  pathwrite = path+"ke.txt"
  griddata = cbinread(path+"all.scalars-bous",-1,"double")
  ndata1 = sizeof(griddata)
  ndata2 =ndata1/8. ; 8 bytes per double
  ndata = doubletoint(ndata2)
  print(ndata)
  print(griddata)
;
; Compute the total number of time increments in the data
;
  icount = 0
  ni = (/griddata(0)/)
  nints = doubletoint(ni)
  if (nints.eq.8) then
    stride=10
  else
    stride=16
  end if
;
;
;
  dims = 0
do while (icount.lt.ndata)
   ni = (/griddata(icount)/)
   nints = doubletoint(ni)
   if(nints.eq.8) then
    stride=10
   else
    stride=16
   end if
   icount = icount + stride
   dims=dims + 1
  end do
;
;
;
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
  pv_qg=new(dim,"double")
  potens_qg=new(dim,"double")
  pv_r0=new(dim,"double")
  potens_r0=new(dim,"double")
  pv_fr0=new(dim,"double")
  potens_fr0=new(dim,"double")
  print(griddata)
;
; defaults
; 
nx = 0
fcor = 0
f_k = 0
;
; Loop over all possible data
;
;
;
; If nints = 8 then start with stride with a factor of 10
;
    ni = (/griddata(0)/)
    nints = doubletoint(ni)
    print(nints)
    time(0) = (/griddata(1)/)
    print(time(0))
;
    if (nints.eq.8) then
       stride=10
    else
       stride=16
    end if
;
;
;
  icount = 0
  stride=0
  do while (icount .le.  dims-1)
;
;   Assign data
;
    print(icount)
    ni = (/griddata(stride)/)
    nints = doubletoint(ni)
    print(nints)
    time(icount) = (/griddata(stride + 1)/)
    print(time(icount))
    ke(icount) = (/griddata(stride + 2)/)
    pe(icount) = (/griddata(stride + 3)/)
    ke_diss(icount) = (/griddata(stride + 4)/)
    pe_diss(icount) = (/griddata(stride + 5)/)
    pv(icount) = (/griddata(stride + 6)/)
    pv_diss(icount) = (/griddata(stride + 7)/)
    potens(icount) = (/griddata(stride + 8)/)
    enstr_diss(icount) = (/griddata(stride + 9)/)
    if (nints.eq.8) then
       pv_qg(icount) = 0.
       potens_qg(icount) = 0.
       pv_r0(icount) = 0.
       potens_r0(icount) = 0.
       pv_fr0(icount) = 0.
       potens_fr0(icount) = 0.
    else
       pv_qg(icount) = (/griddata(stride + 10)/)
       potens_qg(icount) = (/griddata(stride + 11)/)
       pv_r0(icount) = (/griddata(stride + 12)/)
       potens_r0(icount) = (/griddata(stride + 13)/)
       pv_fr0(icount) = (/griddata(stride + 14)/)
       potens_fr0(icount) = (/griddata(stride + 15)/)
    end if
    if (nints.eq.8) then
       stride=stride + 10
    else
       stride=stride + 16
    end if
    icount = icount + 1
  end do
;
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
;    Write out the qg potential vorticity
;     
     toprint(:,0)=time(:)
     toprint(:,1)=pv_qg(:)
     opt = True
     opt@fout = path+"potvor_qg.txt"
     write_matrix(toprint, "2e15.6", opt)
;
;    Write out the qg potential enstrophy 
;     
     toprint(:,0)=time(:)
     toprint(:,1)=potens_qg(:)
     opt = True
     opt@fout = path+"potens_qg.txt"
     write_matrix(toprint, "2e15.6", opt)
;
;    Write out the r0 potential vorticity
;     
     toprint(:,0)=time(:)
     toprint(:,1)=pv_r0(:)
     opt = True
     opt@fout = path+"potvor_r0.txt"
     write_matrix(toprint, "2e15.6", opt)
;
;    Write out the r0 potential enstrophy 
;     
     toprint(:,0)=time(:)
     toprint(:,1)=potens_r0(:)
     opt = True
     opt@fout = path+"potens_r0.txt"
     write_matrix(toprint, "2e15.6", opt)
;
;    Write out the fr0 potential vorticity
;     
     toprint(:,0)=time(:)
     toprint(:,1)=pv_fr0(:)
     opt = True
     opt@fout = path+"potvor_fr0.txt"
     write_matrix(toprint, "2e15.6", opt)
;
;    Write out the fr0 potential enstrophy 
;     
     toprint(:,0)=time(:)
     toprint(:,1)=potens_fr0(:)
     opt = True
     opt@fout = path+"potens_fr0.txt"
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
