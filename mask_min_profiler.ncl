;
; File: mask_min_profiler.ncl
;
; Author: D. Adriaansen
;
; Date: 05 May 2016
;
; Purpose: Look at the precip flag and determine the number of periods longer than a certain
;          time threshold defined by the user and count them for each day. Then mask out the data
;          and write them out to a file.
;
; Notes:
;
;________________________________________________________________________________________

;#################################### User Config ######################################;

; Path to input netCDF files
ncpath = "/d1/dadriaan/paper/data/masked"

; Path to output netCDF files
ncout = "/d1/dadriaan/paper/data/maskedmin"

; What's the minimum number of minutes you'd like a "good" period to be?
if .not.isvar("mingood") then
  mingood = 240 ; 4 hrs
end if

;###########################################################################################;

; NCL Libraries
load "/home/dadriaan/cvs/pat/dadriaan/scripts/ncl/load_ncl_libs.ncl"
load "/home/dadriaan/cvs/pat/dadriaan/scripts/ncl/ncl_usr_funcs.ncl"

; Open all the files
allfiles = systemfunc("ls "+ncpath+"/*2006*.nc")
fin = addfiles(allfiles,"r")
ListSetType(fin,"cat")

; Read the data
bf = fin[0]->bf
pagl = fin[0]->pagl
palt = fin[0]->palt
plon = fin[0]->plon
plat = fin[0]->plat

unix_time = fin[:]->unix_time
precipflag = fin[:]->precipflag
omegpass1 = fin[:]->omegpass1
zonalpass1 = fin[:]->zonalpass1
meridpass1 = fin[:]->meridpass1
rainpass1 = fin[:]->rainpass1

; Determine the number of times and heights
ntimes = dimsizes(unix_time)
numhgts = dimsizes(pagl)

; Print info for user
nhrsdata = dimsizes(unix_time)/60
ndaysdata = dimsizes(unix_time)/1440
print("PROCESSING "+ndaysdata+" DAYS OF DATA ("+nhrsdata+" HOURS)")
print("TOTAL MINUTES = "+dimsizes(unix_time))
print("")

print("=================PROCESSING ALL")
print("")
; Find the indices of all the times identified as precip
badtimes = ind(precipflag.gt.0)
print("FOUND "+dimsizes(badtimes)+" times with precip.")

; Loop over the bad times, find the length of each bad period and store it
badprecip = new((/dimsizes(badtimes)/),integer)
do f=0,dimsizes(badtimes)-2,1
  dnt = badtimes(f+1)-badtimes(f)
  if dnt.gt.1 then
    nmin = dnt%60
    nhrs = dnt/60
    ;print("")
    ;print("LENGTH OF PERIOD = "+nhrs+" HRS "+nmin+" MIN")
    ;print("PER BEG IDX = "+badtimes(f))
    ;print("PER END IDX = "+badtimes(f+1))
    badprecip(f) = dnt
  end if
end do
npers = dimsizes(ind(.not.ismissing(badprecip)))
print("FOUND "+npers+" GOOD PERIODS.")
print("")

; At this point the 50 MHz data have been flagged for precipitation using the 920 data. We now want to
; bring in the idea of a "darwin" day, which runs 02Z to 02Z. But first we want to test some interpolation.
; We need to find periods that were flagged for precip less than some number of minutes (configurable), and
; do linear interpolation between the gaps. We should test 0, 1, 5, and 10 minutes and plot time series of the
; data to make sure it looks "ok" after interpolation.

; Loop over the indices and check how many times each period is. If it's less than some threshold
; then just turn on the precip flag those times too, effectively "filtering" those times.
; The time resolution of the data here is assumed to be 1 minute.
print("===============ELIMINATING MINGOOD")
print("")
do f=0,dimsizes(badtimes)-2,1
  dnt = badtimes(f+1)-badtimes(f)
  ;print("NUM TIMES = "+dnt)
  ;if dnt.gt.1 .and. dnt.le.mingood then
  if dnt.le.mingood then
    precipflag(badtimes(f):badtimes(f+1)) = 1
    ;print("NUM TIMES = "+dnt)
  end if
  ;if dnt.gt.1 then
  ;  print("NUM ANY TIMES = "+dnt)
  ;end if
end do

; Now do it again and print the information about the length of the good periods.
newbad = ind(precipflag.gt.0)
badminprecip = new((/dimsizes(newbad)/),integer)
print("FOUND "+dimsizes(newbad)+" total 'bad' (precip + less than min) times.")
print("")
do f=0,dimsizes(newbad)-2,1
  dnt = newbad(f+1)-newbad(f)
  if dnt.gt.1 then
    ;print("NUM NEW TIMES = "+dnt)
    nmin = dnt%60
    nhrs = dnt/60
    ;print("")
    ;print("LENGTH OF PERIOD = "+nhrs+" HRS "+nmin+" MIN")
    ;print("PER BEG IDX = "+newbad(f))
    ;print("PER END IDX = "+newbad(f+1))
    badminprecip(f) = dnt
  end if
end do
npers2 = dimsizes(ind(.not.ismissing(badminprecip)))
print("FOUND "+npers2+" GOOD PERIODS.")
print("")

;========================
; USE CODE BELOW TO APPLY THE NEW MASK TO THE DATA
;
; Store the final pass after filtering out "noisy" periods of good data less than threshold
rainpass2 = new((/dimsizes(rainpass1)/),float)
rainpass2@_FillValue = bf
omegpass2 = new((/dimsizes(omegpass1)/),float)
omegpass2@_FillValue = bf
zonalpass2 = new((/dimsizes(zonalpass1)/),float)
zonalpass2@_FillValue = bf
meridpass2 = new((/dimsizes(meridpass1)/),float)
meridpass2@_FillValue = bf
do nt=0,ntimes-1,1
  if precipflag(nt).gt.0 then
    rainpass2(nt,:) = bf
    omegpass2(nt,:) = bf
    zonalpass2(nt,:) = bf
    meridpass2(nt,:) = bf
  else
    rainpass2(nt,:) = rainpass1(nt,:)
    omegpass2(nt,:) = omegpass1(nt,:)
    zonalpass2(nt,:) = zonalpass1(nt,:)
    meridpass2(nt,:) = meridpass1(nt,:)
  end if
end do

; Determine the number of output files to write
nof = dimsizes(allfiles)
do f=0,nof-1,1
  
  ; Determine the beg and end to write data for
  beg_ind = f*1440
  end_ind = beg_ind+1439
  ntimes = (end_ind-beg_ind)+1
  
  ; Set the filename
  s1 = str_split(allfiles(f),"/")
  s2 = s1(dimsizes(s1)-1)
  fout = ncout+"/maskmin."+s2
  
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
  filevardef(fo,"unix_time",typeof(unix_time),(/"time"/))
  filevardef(fo,"bf",typeof(bf),(/"single_value"/))
  filevardef(fo,"palt",typeof(palt),(/"single_value"/))
  filevardef(fo,"plat",typeof(plat),(/"single_value"/))
  filevardef(fo,"plon",typeof(plon),(/"single_value"/))
  filevardef(fo,"pagl",typeof(pagl),(/"height"/))
  filevardef(fo,"precipflag",typeof(precipflag),(/"time"/))
  filevardef(fo,"omegpass1",typeof(omegpass1),(/"time","height"/))
  filevardef(fo,"zonalpass1",typeof(zonalpass1),(/"time","height"/))
  filevardef(fo,"meridpass1",typeof(meridpass1),(/"time","height"/))
  filevardef(fo,"rainpass1",typeof(rainpass1),(/"time","height"/))
  filevardef(fo,"omegpass2",typeof(omegpass2),(/"time","height"/))
  filevardef(fo,"zonalpass2",typeof(zonalpass2),(/"time","height"/))
  filevardef(fo,"meridpass2",typeof(meridpass2),(/"time","height"/))
  filevardef(fo,"rainpass2",typeof(rainpass2),(/"time","height"/))

  ; Define output variable attributes
  filevarattdef(fo,"unix_time",unix_time)
  filevarattdef(fo,"bf",bf)
  filevarattdef(fo,"palt",palt)
  filevarattdef(fo,"plat",plat)
  filevarattdef(fo,"plon",plon)
  filevarattdef(fo,"pagl",pagl)
  filevarattdef(fo,"precipflag",precipflag)
  filevarattdef(fo,"omegpass1",omegpass1)
  filevarattdef(fo,"zonalpass1",zonalpass1)
  filevarattdef(fo,"meridpass1",meridpass1)
  filevarattdef(fo,"rainpass1",rainpass1)
  filevarattdef(fo,"omegpass2",omegpass2)
  filevarattdef(fo,"zonalpass2",zonalpass2)
  filevarattdef(fo,"meridpass2",meridpass2)
  filevarattdef(fo,"rainpass2",rainpass2)

  ; Exit define mode
  setfileoption(fo,"DefineMode",False)

  ; Write data to the file
  fo->unix_time = (/unix_time(beg_ind:end_ind)/)
  fo->bf = (/bf/)
  fo->palt = (/palt/)
  fo->plat = (/plat/)
  fo->plon = (/plon/)
  fo->pagl = (/pagl/)
  fo->precipflag = (/precipflag(beg_ind:end_ind)/)
  fo->omegpass1 = (/omegpass1(beg_ind:end_ind,:)/)
  fo->zonalpass1 = (/zonalpass1(beg_ind:end_ind,:)/)
  fo->meridpass1 = (/meridpass1(beg_ind:end_ind,:)/)
  fo->rainpass1 = (/rainpass1(beg_ind:end_ind,:)/)
  fo->omegpass2 = (/omegpass2(beg_ind:end_ind,:)/)
  fo->zonalpass2 = (/zonalpass2(beg_ind:end_ind,:)/)
  fo->meridpass2 = (/meridpass2(beg_ind:end_ind,:)/)
  fo->rainpass2 = (/rainpass2(beg_ind:end_ind,:)/)

end do
