;
; File: mask_50MHz_profiler.ncl
;
; Author: D. Adriaansen
;
; Date: 29 Apr 2016
;
; Purpose: mask the 50 MHz profiler data with 920 MHz information (precip)
;          and save to netCDF files
;
; Notes:
;______________________________________________________________________________

;############################### User Config #####################################;

; Path to input netCDF files
ncpath = "/d1/dadriaan/paper/data/raw"

; Path to output netCDF files
ncout = "/d1/dadriaan/paper/data/masked_offset"

; Which tests would you like to employ to filter te 920 data for precip?
; TEST 1: z > 4 km && valid dbz
t1 = True
t1z = 4000.0

; TEST 2: abs(dopvel) > 1.0 m/s
t2 = False
t2v = 10.0

; TEST 3: dbz > 20.0 dbz
t3 = False
t3r = 20.0

; Which tests would you like to employ to filter the 50 MHz data for clouds?
; TEST4: abs(omega) > 1.0 m/s and Z < 8 km
t4 = False
t4v = 5.0
t4z = 8000.0

; What hour is the beginning of a day? In Darwin, we will use 02Z to 02Z, or 1130-1130 local time.
beghr = 2

; Print configuration for user:
print("")
print("#############################")
print("")
print("DAY BEGINS AT: "+sprinti("%02i",tointeger(beghr))+"Z")
;if daynum.lt.0 then
;  print("PLOTTING DATA FOR ALL DAYS")
;else
;  print("PLOTTING DATA FOR DAY NUMBER: "+daynum)
;end if
if t1 then
  print("PRECIP TEST 1: YES")
else
  print("PRECIP TEST 2: NO")
end if
if t2 then
  print("PRECIP TEST 2: YES")
else
  print("PRECIP TEST 2: NO")
end if
if t3 then
  print("PRECIP TEST 3: YES")
else
  print("PRECIP TEST 3: NO")
end if
print("")
print("############################")
print("")

;###################################################################################;

; NCL Libraries
load "/home/dadriaan/cvs/pat/dadriaan/scripts/ncl/load_ncl_libs.ncl"
load "/home/dadriaan/cvs/pat/dadriaan/scripts/ncl/ncl_usr_funcs.ncl"

; Open the netcdf files
allfiles = systemfunc("ls "+ncpath+"/*2006*.nc")
fin = addfiles(allfiles,"r")
ListSetType(fin,"cat")

; Read the data only needed from a single file
bf = fin[0]->bad_flag
palt = fin[0]->prof_site_altitude
plat = fin[0]->prof_site_latitude
plon = fin[0]->prof_site_longitude
pagl = fin[0]->prof_height_100m_AGL

; Time information
y = fin[:]->prof_time_year
m = fin[:]->prof_time_month
jd = fin[:]->prof_time_dayofyear
d = fin[:]->prof_time_dayofmonth
h = fin[:]->prof_time_hour
mn = fin[:]->prof_time_minute

; Compute some info for the user
nhrsdata = dimsizes(y)/60
ndaysdata = dimsizes(y)/1440
print("PROCESSING "+ndaysdata+" DAYS OF DATA ("+nhrsdata+" HOURS)")
print("TOTAL MINUTES = "+dimsizes(y))
print("")

; Read in 50 MHz data and assign missing data code
omean = fin[:]->prof_dar50_Omega_mean
omean@_FillValue = bf
zmean = fin[:]->prof_dar50_Zonal_mean
zmean@_FillValue = bf
mmean = fin[:]->prof_dar50_Merid_mean
mmean@_FillValue = bf

; Read in 920 MHz data and assign missing data code
rainsnr = fin[:]->prof_dar920_vert_snr
rainsnr@_FillValue = bf
rainvel = fin[:]->prof_dar920_vert_vel
rainvel@_FillValue = bf
rainwid = fin[:]->prof_dar920_vert_wid
rainwid@_FillValue = bf
raindbz = fin[:]->prof_dar920_vert_zdb
raindbz@_FillValue = bf
rainhfs = fin[:]->prof_hydrometeor_fall_speed
rainhfs@_FillValue = bf

; Determine number of times and heights
ntimes = dimsizes(y)
numhgts = dimsizes(pagl)

; Her array to hold times
utimes = new((/ntimes/),integer)
utimes@units = "seconds since 1970-01-01 00:00:00"

; New array for times to ignore
precipflag = new((/ntimes/),integer)
precipflag(:) = 0

; Loop over the time and create a continuous UNIX time variable
do t=0,ntimes-1,1
  utimes(t) = toint(cd_inv_calendar(toint(y(t)),toint(m(t)),toint(d(t)),toint(h(t)),toint(mn(t)),0,"seconds since 1970-1-1 00:00:00",0))
end do

; Set up 2D fields of U,V,W 50 MHz data indicating whether the data are good, bad, or "precip"
; Initialize to 1 (assume all points are good initially, good data = 1)
mask_w = new((/dimsizes(omean)/),integer)
mask_u = new((/dimsizes(zmean)/),integer)
mask_v = new((/dimsizes(mmean)/),integer)
mask_w@_FillValue = -9
mask_u@_FillValue = -9
mask_v@_FillValue = -9
mask_w(:,:) = 1
mask_v(:,:) = 1
mask_u(:,:) = 1
mask_w!0 = "time"
mask_w!1 = "height"
mask_w@height = pagl
mask_u!0 = "time"
mask_u!1 = "height"
mask_u@height = pagl
mask_v!0 = "time"
mask_v!1 = "height"
mask_v@height = pagl

; Set up 2D fields of 920 data indicating whether the data are good, bad, or "precip"
; Initialize to 1 (assume all points are good initially, good data = 1)
mask_dbz = new((/dimsizes(raindbz)/),integer)
mask_dbz@_FillValue = -9
mask_dbz(:,:) = 1
mask_dbz!0 = "time"
mask_dbz!1 = "height"
mask_dbz@height = pagl

; Loop to set the values of the mask fields based on the bad data flag (bad data = 2)
do nz=0,numhgts-1,1
  miss_w = ind(ismissing(omean(:,nz)))
  if .not.ismissing(miss_w(0)) .and. dimsizes(miss_w).gt.0 then
    mask_w(miss_w,nz) = 2
  end if
  miss_u = ind(ismissing(zmean(:,nz)))
  if .not.ismissing(miss_u(0)) .and. dimsizes(miss_u).gt.0 then
    mask_u(miss_u,nz) = 2
  end if
  miss_v = ind(ismissing(mmean(:,nz)))
  if .not.ismissing(miss_v(0)) .and. dimsizes(miss_v).gt.0 then
    mask_v(miss_v,nz) = 2
  end if
  miss_p = ind(ismissing(raindbz(:,nz)))
  if .not.ismissing(miss_p(0)) .and. dimsizes(miss_p).gt.0 then
    mask_dbz(miss_p,nz) = 2
  end if
  delete([/miss_w,miss_u,miss_v,miss_p/])
end do

; Loop over the data and implement the filtering using the 920 data
do nt=0,ntimes-1,1
  ;print("")
  ;print("PROCESSING TIME: "+utimes(nt))
  do nz=0,numhgts-1,1
    ; TEST 1
    if t1 then
      if pagl(nz).gt.t1z .and. .not.ismissing(raindbz(nt,nz)) then
        precipflag(nt) = 1
        ;print("")
        ;print("TEST 1 TRIGGER")
        ;print("pagl(nz) = "+pagl(nz))
        ;print("raindbz(nz) = "+raindbz(nt,nz))
        break
      end if
    end if

    ; TEST 2
    if t2 then
      if .not.ismissing(rainvel(nt,nz)) .and. abs(rainvel(nt,nz)).gt.t2v then
        precipflag(nt) = 1
        ;print("")
        ;print("TEST 2 TRIGGER")
        ;print("rainvel(nz) = "+rainvel(nt,nz))
        break
      end if
    end if

    ; TEST 3
    if t3 then
      if .not.ismissing(raindbz(nt,nz)) .and. raindbz(nt,nz).gt.t3r then
        precipflag(nt) = 1
        ;print("")
        ;print("TEST 3 TRIGGER")
        ;print("raindbz(nz) = "+raindbz(nt,nz))
        break
      end if
    end if

    ; TEST 4
    if t4 then
      if .not.ismissing(omean(nt,nz)) .and. (abs(omean(nt,nz)).gt.t4v .and. pagl(nz).lt.t4z)then
        precipflag(nt) = 1
        ;print("")
        ;print("TEST 4 TRIGGER")
        ;print("abs(omean(nz)) = "+abs(omean(nt,nz)))
        break
      end if
    end if
  end do
end do

; Now use the 1D precipitation flag vector to mask out the 2D data (precip = 3)
do nt=0,ntimes-1,1
  if precipflag(nt).gt.0 then
    mask_w(nt,:) = 3
    mask_v(nt,:) = 3
    mask_u(nt,:) = 3
    mask_dbz(nt,:) = 3
  end if
end do

; Apply the precip flag to the data
rainpass1 = new((/dimsizes(raindbz)/),float)
rainpass1@_FillValue = bf
omegpass1 = new((/dimsizes(omean)/),float)
omegpass1@_FillValue = bf
zonalpass1 = new((/dimsizes(zmean)/),float)
zonalpass1@_FillValue = bf
meridpass1 = new((/dimsizes(mmean)/),float)
meridpass1@_FillValue = bf
do nt=0,ntimes-1,1
  if precipflag(nt).gt.0 then
    rainpass1(nt,:) = bf
    omegpass1(nt,:) = bf
    zonalpass1(nt,:) = bf
    meridpass1(nt,:) = bf
  else
    rainpass1(nt,:) = raindbz(nt,:)
    omegpass1(nt,:) = omean(nt,:)
    zonalpass1(nt,:) = zmean(nt,:)
    meridpass1(nt,:) = mmean(nt,:)
  end if
end do

; What variables do we want to write out to the file?
; bf, palt, plat, plon, pagl
; utimes, precipflag
; omegpass1, meridpass1, zonalpass1, rainpass1
; mask_u, mask_w, mask_v, mask_dbz
; omean, zmean, mmean, raindbz
; To keep simplicity, we will write out data in the same way it came in- days at a time starting at 00Z

; Determine number of output files to write
nof = dimsizes(allfiles)
do f=0,nof-1,1

  ; Determine the beg and end to write data for
  beg_ind = f*1440
  end_ind = beg_ind+1439
  ntimes = (end_ind-beg_ind)+1
  ;print("NUM TIMES = "+ntimes)

  ; Set the filename
  s1 = str_split(allfiles(f),"/")
  s2 = s1(dimsizes(s1)-1)
  fout = ncout+"/mask."+s2

  ; Print for user
  print("WRITING FILE: "+fout)
  print("")

  ; Prepare for writing
  system("/bin/rm -f "+fout)           ; Remove any existing file with the same output name
  fo = addfile(fout,"c")               ; Create the netcdf file
  setfileoption(fo,"DefineMode",True)  ; Enter define mode
  fAtt = True                          ; Set file attributes to True
  fAtt@title = "Precipitation masked 50 MHz and 920 MHz data"        ; Set the title of the file
  fileattdef(fo,fAtt)                  ; Add the attributes to the output file

  ; Output file dimensions
  dimNames = (/"time","height","single_value"/)
  dimSizes = (/ntimes,numhgts,1/)
  dimUnlim = (/False,False,False/)

  ; Define output file dimensions
  filedimdef(fo,dimNames,dimSizes,dimUnlim)
  
  ; Define output variables
  filevardef(fo,"unix_time",typeof(utimes),(/"time"/))
  filevardef(fo,"bf",typeof(bf),(/"single_value"/))
  filevardef(fo,"palt",typeof(palt),(/"single_value"/))
  filevardef(fo,"plat",typeof(plat),(/"single_value"/))
  filevardef(fo,"plon",typeof(plon),(/"single_value"/))
  filevardef(fo,"pagl",typeof(pagl),(/"height"/))
  filevardef(fo,"omean",typeof(omean),(/"time","height"/))
  filevardef(fo,"zmean",typeof(zmean),(/"time","height"/))
  filevardef(fo,"mmean",typeof(mmean),(/"time","height"/))
  filevardef(fo,"raindbz",typeof(raindbz),(/"time","height"/))
  filevardef(fo,"mask_w",typeof(mask_w),(/"time","height"/))
  filevardef(fo,"mask_u",typeof(mask_u),(/"time","height"/))
  filevardef(fo,"mask_v",typeof(mask_v),(/"time","height"/))
  filevardef(fo,"mask_dbz",typeof(mask_dbz),(/"time","height"/))
  filevardef(fo,"precipflag",typeof(precipflag),(/"time"/))
  filevardef(fo,"omegpass1",typeof(omegpass1),(/"time","height"/))
  filevardef(fo,"zonalpass1",typeof(zonalpass1),(/"time","height"/))
  filevardef(fo,"meridpass1",typeof(meridpass1),(/"time","height"/))
  filevardef(fo,"rainpass1",typeof(rainpass1),(/"time","height"/))
  
  ; Define output variable attributes
  filevarattdef(fo,"unix_time",utimes)
  filevarattdef(fo,"bf",bf)
  filevarattdef(fo,"palt",palt)
  filevarattdef(fo,"plat",plat)
  filevarattdef(fo,"plon",plon)
  filevarattdef(fo,"pagl",pagl)
  filevarattdef(fo,"omean",omean)
  filevarattdef(fo,"zmean",zmean)
  filevarattdef(fo,"mmean",mmean)
  filevarattdef(fo,"raindbz",raindbz)
  filevarattdef(fo,"mask_w",mask_w)
  filevarattdef(fo,"mask_u",mask_u)
  filevarattdef(fo,"mask_v",mask_v)
  filevarattdef(fo,"mask_dbz",mask_dbz)
  filevarattdef(fo,"omegpass1",omegpass1)
  filevarattdef(fo,"zonalpass1",zonalpass1)
  filevarattdef(fo,"meridpass1",meridpass1)
  filevarattdef(fo,"rainpass1",rainpass1)
  
  ; Exit define mode
  setfileoption(fo,"DefineMode",False)
  
  ; Write data to the file
  fo->unix_time = (/utimes(beg_ind:end_ind)/)
  fo->bf = (/bf/)
  fo->palt = (/palt/)
  fo->plat = (/plat/)
  fo->plon = (/plon/)
  fo->pagl = (/pagl/)
  fo->omean = (/omean(beg_ind:end_ind,:)/)
  fo->zmean = (/zmean(beg_ind:end_ind,:)/)
  fo->mmean = (/mmean(beg_ind:end_ind,:)/)
  fo->raindbz = (/raindbz(beg_ind:end_ind,:)/)
  fo->mask_w = (/mask_w(beg_ind:end_ind,:)/)
  fo->mask_u = (/mask_u(beg_ind:end_ind,:)/)
  fo->mask_v = (/mask_v(beg_ind:end_ind,:)/)
  fo->mask_dbz = (/mask_dbz(beg_ind:end_ind,:)/)
  fo->omegpass1 = (/omegpass1(beg_ind:end_ind,:)/)
  fo->zonalpass1 = (/zonalpass1(beg_ind:end_ind,:)/)
  fo->meridpass1 = (/meridpass1(beg_ind:end_ind,:)/)
  fo->rainpass1 = (/rainpass1(beg_ind:end_ind,:)/)

end do
