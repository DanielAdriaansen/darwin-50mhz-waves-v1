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
ncout = "/d1/dadriaan/paper/data/maskedminbad"

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
mask_w = fin[:]->mask_w
mask_u = fin[:]->mask_u
mask_v = fin[:]->mask_v
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

; New matrices to hold the final filtered data
rainpass2 = new((/dimsizes(rainpass1)/),float)
rainpass2@_FillValue = bf
omegpass2 = new((/dimsizes(omegpass1)/),float)
omegpass2@_FillValue = bf
zonalpass2 = new((/dimsizes(zonalpass1)/),float)
zonalpass2@_FillValue = bf
meridpass2 = new((/dimsizes(meridpass1)/),float)
meridpass2@_FillValue = bf

print("=================PROCESSING ALL WINDS")
print("")
; This needs to operate in 3D- do it level by level but also component by component since
; sometimes w could be good but u or v might be bad. To find bad data, use the "mask_x" matrices
; which contain flags for the data quality (1 = good, 2 = bad, 3 = precip)
;do z=29,29,1
do z=0,numhgts-1,1

  ; Print current height for user
  print("z = "+pagl(z))
  
  ; Find the indices of all the times identified as precip or bad data at this height
  badw = ind(mask_w(:,z).gt.1) ; PRECIP + BAD
  badu = ind(mask_u(:,z).gt.1) ; PRECIP + BAD
  badv = ind(mask_v(:,z).gt.1) ; PRECIP + BAD
  ;badw = ind(mask_w(:,z).eq.3) ; PRECIP ONLY
  ;badu = ind(mask_u(:,z).eq.3) ; PRECIP ONLY
  ;badv = ind(mask_v(:,z).eq.3) ; PRECIP ONLY
  print("FOUND "+dimsizes(badw)+" times that don't work (precip or bad) for W.")
  print("PERCENT TIMES IGNORED FOR W = "+(tofloat(dimsizes(badw))/tofloat(ntimes))*100.0)
  print("FOUND "+dimsizes(badu)+" times that don't work (precip or bad) for U.")
  print("PERCENT TIMES IGNORED FOR U = "+(tofloat(dimsizes(badu))/tofloat(ntimes))*100.0)
  print("FOUND "+dimsizes(badv)+" times that don't work (precip or bad) for V.")
  print("PERCENT TIMES IGNORED FOR V = "+(tofloat(dimsizes(badv))/tofloat(ntimes))*100.0)

  ; Loop over the bad times (indexes), find the length of each period and store it. Since the data
  ; are spaced every 1 minute, using the indexes is the same as using minutes. Taking the difference
  ; between adjacent indexes will identify the length (minutes) of the "good" periods of data. For example, if
  ; the array of indexes jumps from 140 to 200, we know there is a 58 minute period of "good" data because
  ; we indexed looking for bad data. The "good" data period would be from index 141 to 199.
  print("=============FINDING GOOD PERIODS.")
  badwdata = new((/dimsizes(badw)/),integer)
  badudata = new((/dimsizes(badu)/),integer)
  badvdata = new((/dimsizes(badv)/),integer)
  do f=0,dimsizes(badw)-2,1
    dnt = badw(f+1)-badw(f)
    if dnt.gt.1 then
      gbeg = badw(f)+1
      gend = badw(f+1)-1
      ;print("GBEG = "+gbeg)
      ;print("GEND = "+gend)
      gdiff = gend-gbeg
      if gdiff.eq.0.0 then
        ;print("FOUND 0")
        ;print("badf = "+badw(f))
        ;print("badf+1 = "+badw(f+1))
        ;print("gbeg = "+gbeg)
        ;print("gend = "+gend)
        nmin = 1
        nhrs = 0
      else
        ;print("FOUND >0")
        ;print("badf = "+badw(f))
        ;print("badf+1 = "+badw(f+1))
        ;print("gbeg = "+gbeg)
        ;print("gend = "+gend)
        nmin = gdiff%60
        nhrs = gdiff/60
      end if
      ;;nmin = dnt%60
      ;;nhrs = dnt/60
      ;print("")
      ;print("LENGTH OF PERIOD = "+nhrs+" HRS "+nmin+" MIN")
      ;print("PER BEG IDX = "+badw(f))
      ;print("PER END IDX = "+badw(f+1))
      ;print("PER BEG IDX = "+gbeg)
      ;print("PER END IDX = "+gend)
      ;badwdata(f) = dnt
      badwdata(f) = nmin
    end if
  end do
  do f=0,dimsizes(badu)-2,1
    dnt = badu(f+1)-badu(f)
    if dnt.gt.1 then
      gbeg = badu(f)+1
      gend = badu(f+1)-1
      gdiff = gend-gbeg
      if gdiff.eq.0.0 then
        nmin = 1
        nhrs = 0
      else
        nmin = gdiff%60
        nhrs = gdiff/60
      end if
      ;badudata(f) = dnt
      badudata(f) = nmin
    end if
  end do
  do f=0,dimsizes(badv)-2,1
    dnt = badv(f+1)-badv(f)
    if dnt.gt.1 then
      gbeg = badv(f)+1
      gend = badv(f+1)-1
      gdiff = gend-gbeg
      if gdiff.eq.0.0 then
        nmin = 1
        nhrs = 0
      else
        nmin = gdiff%60
        nhrs = gdiff/60
      end if
      ;badvdata(f) = dnt
      badvdata(f) = nmin
    end if
  end do
  ; Print the number of valid good periods. This is simply anywhere the array baddata has valid data, since
  ; baddata contains the length of each good period in minutes.
  npers = dimsizes(ind(.not.ismissing(badwdata)))
  print("FOUND "+npers+" GOOD W PERIODS BEFORE FILTER.")
  npers = dimsizes(ind(.not.ismissing(badudata)))
  print("FOUND "+npers+" GOOD U PERIODS BEFORE FILTER.")
  npers = dimsizes(ind(.not.ismissing(badvdata)))
  print("FOUND "+npers+" GOOD V PERIODS BEFORE FILTER.")

  ;/// JUNK COMMENT REMOVE EVENTUALLY
  ; At this point the 50 MHz data have been flagged for precipitation using the 920 data. We now want to
  ; bring in the idea of a "darwin" day, which runs 02Z to 02Z. But first we want to test some interpolation.
  ; We need to find periods that were flagged for precip less than some number of minutes (configurable), and
  ; do linear interpolation between the gaps. We should test 0, 1, 5, and 10 minutes and plot time series of the
  ; data to make sure it looks "ok" after interpolation.
  ;///

  ; Now do the same looping again, but "turn off" any good periods that are less than the threshold "mingood" specified
  ; by the user at the top of this script. The time resolution of the data here is assumed to be 1 minute.
  print("===============ELIMINATING MINGOOD")
  do f=0,dimsizes(badw)-2,1
    dnt = badw(f+1)-badw(f)
    gbeg = badw(f)+1
    gend = badw(f+1)-1
    glen = gend-gbeg
    ;print("NUM TIMES = "+dnt)
    if dnt.gt.1 .and. glen.lt.mingood then
    ;if glen.le.mingood then
      ;precipflag(badtimes(f):badtimes(f+1)) = 1
      ;mask_w(badw(f):badw(f+1),z) = 4 ; ASSIGN SPECIAL VALUE OF TOO SMALL
      mask_w(gbeg:gend,z) = 4 ; ASSIGN SPECIAL VALUE OF TOO SMALL
      ;print("NUM TIMES = "+dnt)
    end if
    ;if dnt.gt.1 then
    ;  print("NUM ANY TIMES = "+dnt)
    ;end if
  end do
  do f=0,dimsizes(badu)-2,1
    dnt = badu(f+1)-badu(f)
    gbeg = badu(f)+1
    gend = badu(f+1)-1
    glen = gend-gbeg
    ;if glen.le.mingood then
    if dnt.gt.1 .and. glen.lt.mingood then
      ;mask_u(badu(f):badu(f+1),z) = 4 ; ASSIGN SPECIAL VALUE OF TOO SMALL
      mask_u(gbeg:gend,z) = 4 ; ASSIGN SPECIAL VALUE OF TOO SMALL
    end if
  end do
  do f=0,dimsizes(badv)-2,1
    dnt = badv(f+1)-badv(f)
    gbeg = badv(f)+1
    gend = badv(f+1)-1
    glen = gend-gbeg
    ;if glen.le.mingood then
    if dnt.gt.1 .and. glen.lt.mingood then
      ;mask_v(badv(f):badv(f+1),z) = 4 ; ASSIGN SPECIAL VALUE OF TOO SMALL
      mask_v(gbeg:gend,z) = 4 ; ASSIGN SPECIAL VALUE OF TOO SMALL
    end if
  end do

  ; Find the indices of all the times identified as precip, bad, or too short at this height
  newbadw = ind(mask_w(:,z).gt.1) ; PRECIP + BAD
  ;newbadw = ind(mask_w(:,z).gt.2) ; PRECIP ONLY
  print("FOUND "+dimsizes(newbadw)+" times that don't work (precip, bad, or too short) for W.")
  print("PERCENT W TIMES IGNORED = "+(tofloat(dimsizes(newbadw))/tofloat(ntimes))*100.0)
  newbadu = ind(mask_u(:,z).gt.1) ; PRECIP + BAD
  ;newbadu = ind(mask_u(:,z).gt.2) ; PRECIP ONLY
  print("FOUND "+dimsizes(newbadu)+" times that don't work (precip, bad, or too short) for U.")
  print("PERCENT U TIMES IGNORED = "+(tofloat(dimsizes(newbadu))/tofloat(ntimes))*100.0)
  newbadv = ind(mask_v(:,z).gt.1) ; PRECIP + BAD
  ;newbadv = ind(mask_v(:,z).gt.2) ; PRECIP ONLY
  print("FOUND "+dimsizes(newbadv)+" times that don't work (precip, bad, or too short) for V.")
  print("PERCENT V TIMES IGNORED = "+(tofloat(dimsizes(newbadv))/tofloat(ntimes))*100.0)

  ; Finally, repeat the looping a third time collecting information about the length of the "good" periods after
  ; filtering using the minimum threshold specified by the user. This is the total number of periods at the current height
  ; for all times that do not have precip or bad data and meet the mingood threshold.
  print("===============FINDING FINAL PERIODS.")
  badwmin = new((/dimsizes(newbadw)/),integer)
  badumin = new((/dimsizes(newbadu)/),integer)
  badvmin = new((/dimsizes(newbadv)/),integer)
  do f=0,dimsizes(newbadw)-2,1
    dnt = newbadw(f+1)-newbadw(f)
    if dnt.gt.1 then
      gbeg = newbadw(f)+1
      gend = newbadw(f+1)-1
      gdiff = gend-gbeg
      if gdiff.eq.0.0 then
        print("FOUND =0")
        print("badf = "+newbadw(f))
        print("badf+1 = "+newbadw(f+1))
        print("gbeg = "+gbeg)
        print("gend = "+gend)
        nmin = 1
        nhrs = 0
      else
        print("FOUND >0")
        print("badf = "+newbadw(f))
        print("badf+1 = "+newbadw(f+1))
        print("gbeg = "+gbeg)
        print("gend = "+gend)
        nmin = gdiff%60
        nhrs = gdiff/60
      end if
      ;print("NUM NEW TIMES = "+dnt)
      ;nmin = dnt%60
      ;nhrs = dnt/60
      print("")
      print("LENGTH OF PERIOD = "+nhrs+" HRS "+nmin+" MIN")
      ;print("PER BEG IDX = "+newbadw(f))
      ;print("PER END IDX = "+newbadw(f+1))
      print("PER BEG IDX = "+gbeg)
      print("PER END IDX = "+gend)
      ;badwmin(f) = dnt
      badwmin(f) = nmin
    end if
  end do
  do f=0,dimsizes(newbadu)-2,1
    dnt = newbadu(f+1)-newbadu(f)
    if dnt.gt.1 then
      gbeg = newbadu(f)+1
      gend = newbadu(f+1)-1
      gdiff = gend-gbeg
      if gdiff.eq.0.0 then
        nmin = 1
        nhrs = 0
      else
        nmin = gdiff%60
        nhrs = gdiff/60
      end if
      ;badumin(f) = dnt
      badumin(f) = nmin
    end if
  end do
  do f=0,dimsizes(newbadv)-2,1
    dnt = newbadv(f+1)-newbadv(f)
    if dnt.gt.1 then
      gbeg = newbadv(f)+1
      gend = newbadv(f+1)-1
      gdiff = gend-gbeg
      if gdiff.eq.0.0 then
        nmin = 1
        nhrs = 0
      else
        nmin = gdiff%60
        nhrs = gdiff/60
      end if
      ;badvmin(f) = dnt
      badvmin(f) = nmin
    end if
  end do
  
  ; Print the number of valid good periods. This is simply anywhere the array badmin has valid data, since
  ; badmin contains the length of each good period in minutes.
  npers2 = dimsizes(ind(.not.ismissing(badwmin)))
  print("FOUND "+npers2+" GOOD W PERIODS AFTER FILTER.")
  npers2 = dimsizes(ind(.not.ismissing(badumin)))
  print("FOUND "+npers2+" GOOD U PERIODS AFTER FILTER.")
  npers2 = dimsizes(ind(.not.ismissing(badvmin)))
  print("FOUND "+npers2+" GOOD V PERIODS AFTER FILTER.")

  ; Mask out the data
  print("==============MASKING DATA")
  print("")
  rainpass2(:,z) = where(precipflag.gt.0,bf,rainpass1(:,z)) ; Since we only use 920 for precip, just mask it for precip
  ;omegpass2(:,z) = where(mask_w(:,z).gt.2,bf,omegpass1(:,z)) ; Mask based on precip only
  omegpass2(:,z) = where(mask_w(:,z).gt.1,bf,omegpass1(:,z)) ; Mask based on precip+bad
  ;zonalpass2(:,z) = where(mask_u(:,z).gt.2,bf,zonalpass1(:,z)) ; Mask based on precip only
  zonalpass2(:,z) = where(mask_u(:,z).gt.1,bf,zonalpass1(:,z)) ; Mask based on precip+bad
  ;meridpass2(:,z) = where(mask_v(:,z).gt.2,bf,meridpass1(:,z)) ; Mask based on precip only
  meridpass2(:,z) = where(mask_v(:,z).gt.1,bf,meridpass1(:,z)) ; Mask based on precip+bad

  ; Clean up before next height
  delete([/badw,badu,badv,npers,newbadw,newbadu,newbadv,npers2,badwdata,badudata,badvdata,badwmin,badumin,badvmin/])
end do

;========================
; USE CODE BELOW TO APPLY THE NEW MASK TO THE DATA
;
; Store the final pass after filtering out "noisy" periods of good data less than threshold
;do nt=0,ntimes-1,1
;  if precipflag(nt).gt.0 then
;    rainpass2(nt,:) = bf
;    omegpass2(nt,:) = bf
;    zonalpass2(nt,:) = bf
;    meridpass2(nt,:) = bf
;  else
;    rainpass2(nt,:) = rainpass1(nt,:)
;    omegpass2(nt,:) = omegpass1(nt,:)
;    zonalpass2(nt,:) = zonalpass1(nt,:)
;    meridpass2(nt,:) = meridpass1(nt,:)
;  end if
;end do

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
  filevardef(fo,"mask_w",typeof(mask_w),(/"time","height"/))
  filevardef(fo,"mask_u",typeof(mask_u),(/"time","height"/))
  filevardef(fo,"mask_v",typeof(mask_v),(/"time","height"/))

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
  filevarattdef(fo,"mask_w",mask_w)
  filevarattdef(fo,"mask_u",mask_u)
  filevarattdef(fo,"mask_v",mask_v)

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
  fo->mask_w = (/mask_w(beg_ind:end_ind,:)/)
  fo->mask_u = (/mask_u(beg_ind:end_ind,:)/)
  fo->mask_v = (/mask_v(beg_ind:end_ind,:)/)

end do
