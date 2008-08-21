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
epsil_f = .0253303064
;epsil_f = 0.517275
;timescale = (epsil_f*(2.*3.141592*f_k)^2)^(.33333)
;energyscale = (epsil_f/(2.*3.141592*f_k))^(-2./3.)
timescale=1.
energyscale=1.
;timescale = (epsil_f*(f_k)^2)^(.33333)
;energyscale = (epsil_f/(f_k))^(-2./3.)
print(timescale)
print(energyscale)
;
;**************************************************
; read in binary data
;**************************************************
  setfileoption("bin","ReadByteOrder","LittleEndian")
;  path  = "../../Projects/KH/Boussinesq/128v2/n1f100_slow/"
;  path  = "../../Projects/KH/Boussinesq/256/n1f1000/"
;  path  = "../../Projects/Boussinesq/512/n1f500/tar/"
;  path  = "../../Projects/Boussinesq/256/kf3/rop3fr1/Data/"
;  path  = "../../Projects/Boussinesq/512/kf3/rop2fr1/Data/"
  path  = "../../Projects/Boussinesq/256/kf3/ro1fr1/Data/"
  pathwrite = path+"keh_slow.txt"
  griddata = cbinread(path+"all.scalars-slow",-1,"double")
  ndata1 = sizeof(griddata)
  ndata2 =ndata1/8. ; 8 bytes per double
  ndata = doubletoint(ndata2)
;  print(ndata)
;  print(griddata)
;
; Compute the total number of time increments in the data
;
  icount = 0
  ni = (/griddata(0)/)
  nints = doubletoint(ni)
  if (nints.eq.8) then
    stride=10
  end if
  if (nints.eq.5) then
    stride=7
  end if
  if (nints.eq.6) then
    stride=8
  end if
  if (nints.eq.14) then
    stride=16
  end if

  print(nints)
;
;
;
  dims = 0
do while (icount.lt.ndata)
   ni = (/griddata(icount)/)
   nints = doubletoint(ni)
   if(nints.eq.8) then
    stride=10
   end if
   if(nints.eq.14) then
    stride=16
   end if
   if(nints.eq.5) then
    stride=7
   end if
   if(nints.eq.6) then
    stride=8
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
  keh_slow = new(dim,"double")
  pe_slow = new(dim,"double")
  kev_slow = new(dim,"double")
  potens_slow = new(dim,"double")
  vor_slow = new(dim,"double")
  vert_en_slow = new(dim,"double")
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
    end if
    if (nints.eq.5) then
       stride = 7
    end if
    if (nints.eq.6) then
       stride = 8
    end if
     if (nints.eq.14) then
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
    keh_slow(icount) = (/griddata(stride + 2)/)
    kev_slow(icount) = (/griddata(stride + 3)/)
    pe_slow(icount) = (/griddata(stride + 4)/)
    potens_slow(icount) = (/griddata(stride + 5)/)
    vor_slow(icount) = (/griddata(stride + 6)/)
    vert_en_slow(icount) = (/griddata(stride + 7)/)
    print(vert_en_slow(icount))
    if (nints.eq.8) then
       stride=stride + 10
    end if
    if (nints.eq.14) then
       stride = stride + 16
    end if
    if (nints .eq. 5) then
       stride = stride + 7
    end if
    if (nints .eq. 6) then
       stride = stride + 8
    end if
    icount = icount + 1
  end do
;
  time = time*timescale
  print(time)
  pe_slow = pe_slow*energyscale
  keh_slow = keh_slow*energyscale
  kev_slow = kev_slow*energyscale
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
     toprint(:,1)=keh_slow(:)
     opt = True
     opt@fout = path+"keh_slow.txt"
     write_matrix(toprint, "2e25.12", opt)
;
;    Write out the potential energy
;
     toprint(:,0)=time(:)
     toprint(:,1)=pe_slow(:)
     opt = True
     opt@fout = path+"pe_slow.txt"
     write_matrix(toprint, "2e25.12", opt)
;
;    Write out the total energy
;
     toprint(:,0)=time(:)
     toprint(:,1)=vert_en_slow(:) + keh_slow(:)
     opt = True
     opt@fout = path+"tote_slow.txt"
     write_matrix(toprint, "2e25.12", opt)
;
;    Write out the vertical kinetic energy
;
     toprint(:,0)=time(:)
     toprint(:,1)=kev_slow(:)
     opt = True
     opt@fout = path+"kev_slow.txt"
     write_matrix(toprint, "2e25.12", opt)
;
;    Write out the vertical vorticity
;
     toprint(:,0)=time(:)
     toprint(:,1)=vor_slow(:)
     opt = True
     opt@fout = path+"vor_slow.txt"
     write_matrix(toprint, "2e25.12", opt)
;
;    Write out the potential enstrophy
;     
     toprint(:,0)=time(:)
     toprint(:,1)=potens_slow(:)
     opt = True
     opt@fout = path+"potens_slow.txt"
     write_matrix(toprint, "2e25.12", opt)
;
;    Write out the vertical total energy
;     
     toprint(:,0)=time(:)
     toprint(:,1)=vert_en_slow(:)
     opt = True
     opt@fout = path+"vert_en_slow.txt"
     write_matrix(toprint, "2e25.12", opt)
;
; now plot is ke, pe, tote versus time
;
  y2dims = (/3,dims/)
  print(y2dims)
  y2 = new(y2dims,"double")
  do i = 0,dims-1
   y2(0,i) = keh_slow(i)
   y2(1,i) = kev_slow(i)
   y2(2,i) = keh_slow(i) + kev_slow(i)
  end do
;
  res = True
;  res@xyExplicitLegendLabels = (/"tote","keh_slow","kev_slow"/)
;  res@xyExplicitLabels = (/"A","B","C"/)
  res@tiMainString           = "toteke=blue, keh=red, kev=green"
  res@tiXAxisString          = "time"
  res@xyLineColors = (/2,3,4/) ; define line colors
  res@xyLineThicknesses = (/1.,2.,5./) ; define line thickness, 1.0 is deflt
  res@xyDashPattern = (/0.,0.,0./)
  plot = gsn_xy(wks,time,y2,res)

  delete(plot)

end
