;
; File: plot_period_length.ncl
;
; Author: D. Adriaansen
;
; Date: 05 May 2016
;
; Purpose: Look at the precip flag and determine the number of periods longer than a certain
;          time threshold defined by the user and count them for each day.
;
; Notes:
;
;________________________________________________________________________________________

;#################################### User Config ######################################;

; Path to input netCDF files
ncpath = "/d1/dadriaan/paper/data/maskedmin"

; Path to image output
imgout = "/d1/dadriaan/paper/images/"

; Monsoon days
mbeg = 13 ; Note- actually day 14, but day 0 = day 1 on zpanel plot
mdays = 20
;mend = 32 ; Note- actually day 33, but day 0 = day 1 on zpanel plot

; Break days
bbeg = 36 ; Note- actually day 37, but day 0 = day 1 on zpanel plot
bdays = 23
;bend = 58 ; Note- actually day 59, but day 0 = day 1 on zpanel plot

; What hour is the beginning of a day? In Darwin, we will use 02Z to 02Z, or 1130-1130 local time.
beghr = 2

; What day number do we want to look at? Use any negative number to inicate we want all days
if .not.isvar("daynum") then
  daynum = -1
end if

; What's the minimum number of minutes you'd like a "good" period to be?
if .not.isvar("mingood") then
  mingood = 240 ; 4 hrs
end if

; Print configuration for user:
print("")
print("#############################")
print("")
print("DAY BEGINS AT: "+sprinti("%02i",tointeger(beghr))+"Z")
if daynum.lt.0 then
  print("PLOTTING DATA FOR ALL DAYS")
else
  print("PLOTTING DATA FOR DAY NUMBER: "+daynum)
end if
print("")
print("############################")
print("")

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
unix_time = fin[:]->unix_time
precipflag = fin[:]->precipflag

; Determine the number of times and heights
ntimes = dimsizes(unix_time)
numhgts = dimsizes(pagl)

; Print info for user
nhrsdata = dimsizes(unix_time)/60
ndaysdata = dimsizes(unix_time)/1440
print("PROCESSING "+ndaysdata+" DAYS OF DATA ("+nhrsdata+" HOURS)")
print("TOTAL MINUTES = "+dimsizes(unix_time))
print("")

; Determine the index of the beginning and the end of the monsoon and break periods
; Need to set the beginning index for the break period
begmonsoon = (beghr*60)+(1440*mbeg)
print("begmonsoon = "+begmonsoon)
print("mbegunix = "+unix_time(begmonsoon))
endmonsoon = begmonsoon+(1440*(mdays-1))
print("endmonsoon = "+endmonsoon)
print("mendunix = "+unix_time(endmonsoon))
ndayssoon = ((endmonsoon-begmonsoon)/1440)
print("ndayssoon = "+ndayssoon)
print("")
begbreak = (beghr*60)+(1440*bbeg)
print("begbreak = "+begbreak)
print("bbegunix = "+unix_time(begbreak))
endbreak = begbreak+(1440*(bdays-1))
print("endbreak = "+endbreak)
print("bendunix = "+unix_time(endbreak))
ndaysbreak = ((endbreak-begbreak)/1440)
print("ndaysbreak = "+ndaysbreak)
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
wks = gsn_open_wks("png",imgout+"histall_"+sprinti("%.3i",toint(mingood)))
hres = True
hres@gsnHistogramBinIntervals = (/0,60,120,180,240,300,360,420,480,540,600,660,720,780,840,900,960,1020,1080,1140,1200,1260,1320,1380,1440,1500,1560,1620,1680,1740,1800/)
hres@tmXBLabelAngleF = 45.0
hres@tmXBLabelFontHeightF = 0.012
hres@tiXAxisString = "Minutes"
hres@tiYAxisString = "Count"
hres@gsnHistogramMinMaxBinsOn = True
hres@gsnDraw = False
hres@gsnFrame = False
bins = gsn_histogram(wks,badprecip,hres)
draw(bins)
frame(wks)

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
wks2 = gsn_open_wks("png",imgout+"histmin_"+sprinti("%.3i",toint(mingood)))
bins2 = gsn_histogram(wks2,badminprecip,hres)
draw(bins2)
frame(wks2)

; At this point, precipflag has been masked by the mingood limits, and the monsoon and break periods can
; be extracted. Once those chunks are isloated, we need to run the mingood code again to turn off periods
; at the beginning and/or end of the monsoon/break that are now less than the mingood.
monflag = precipflag(begmonsoon:endmonsoon)
brkflag = precipflag(begbreak:endbreak)

; Determine the length and number of good periods in the monsoon
monbad = ind(monflag.gt.0)
badmonprecip = new((/dimsizes(monbad)/),integer)
print("===================PROCESSING MONSOON DAYS")
print("")
print("FOUND "+dimsizes(monbad)+ " total 'bad' (precip + less than min) times.")
print("")
do f=0,dimsizes(monbad)-2,1
  dnt = monbad(f+1)-monbad(f)
  if dnt.gt.1 then
    nmin = dnt%60
    nhrs = dnt/60
    ;print("")
    ;print("LENGTH OF PERIOD = "+nhrs+" HRS "+nmin+" MIN")
    ;print("PER BEG IDX = "+monbad(f))
    ;print("PER END IDX = "+monbad(f+1))
    badmonprecip(f) = dnt
  end if
end do
npersmon = dimsizes(ind(.not.ismissing(badmonprecip)))
print("FOUND "+npersmon+" PERIODS BEFORE CHECKING MONSOON BEG/END.")
print("")

print("ELIMINATING MINGOOD (MONSOON)")
print("")
do f=0,dimsizes(monbad)-2,1
  dnt = monbad(f+1)-monbad(f)
  ;print("NUM TIMES = "+dnt)
  ;if dnt.gt.1 .and. dnt.le.mingood then
  if dnt.le.mingood then
    monflag(monbad(f):monbad(f+1)) = 1
    ;print("NUM TIMES = "+dnt)
  end if
  ;if dnt.gt.1 then
  ;  print("NUM ANY TIMES = "+dnt)
  ;end if
end do

print("FOUND "+dimsizes(monbad)+ " total 'bad' (precip + less than min) times.")
print("")
do f=0,dimsizes(monbad)-2,1
  dnt = monbad(f+1)-monbad(f)
  if dnt.gt.1 then
    nmin = dnt%60
    nhrs = dnt/60
    ;print("")
    ;print("LENGTH OF PERIOD = "+nhrs+" HRS "+nmin+" MIN")
    ;print("PER BEG IDX = "+monbad(f))
    ;print("PER END IDX = "+monbad(f+1))
    badmonprecip(f) = dnt
  end if
end do
npersmon = dimsizes(ind(.not.ismissing(badmonprecip)))
print("FOUND "+npersmon+" GOOD PERIODS.")
print("")
wks3 = gsn_open_wks("png",imgout+"histmon_"+sprinti("%.3i",toint(mingood)))
bins3 = gsn_histogram(wks3,badmonprecip,hres)
draw(bins3)
frame(wks3)

; Determine the length and number of good periods in the break
brkbad = ind(brkflag.gt.0)
badbrkprecip = new((/dimsizes(brkbad)/),integer)
print("=================PROCESSING BREAK DAYS")
print("")
print("FOUND "+dimsizes(brkbad)+" total 'bad' (precip + less than min) times.")
print("")
do f=0,dimsizes(brkbad)-2,1
  dnt = brkbad(f+1)-brkbad(f)
  if dnt.gt.1 then
    nmin = dnt%60
    nhrs = dnt/60
    ;print("")
    ;print("LENGTH OF PERIOD = "+nhrs+" HRS "+nmin+" MIN")
    ;print("PER BEG IDX = "+brkbad(f))
    ;print("PER END IDX = "+brkbad(f+1))
    badbrkprecip(f) = dnt
  end if
end do
npersbrk = dimsizes(ind(.not.ismissing(badbrkprecip)))
print("FOUND "+npersbrk+" PERIODS BEFORE CHECKING BREAK BEG/END.")
print("")

print("ELIMINATING MINGOOD (BREAK)")
print("")
do f=0,dimsizes(brkbad)-2,1
  dnt = brkbad(f+1)-brkbad(f)
  ;print("NUM TIMES = "+dnt)
  ;if dnt.gt.1 .and. dnt.le.mingood then
  if dnt.le.mingood then
    brkflag(brkbad(f):brkbad(f+1)) = 1
    ;print("NUM TIMES = "+dnt)
  end if
  ;if dnt.gt.1 then
  ;  print("NUM ANY TIMES = "+dnt)
  ;end if
end do

; Determine the length and number of good periods in the break
brkbad = ind(brkflag.gt.0)
badbrkprecip = new((/dimsizes(brkbad)/),integer)
print("PROCESSING BREAK DAYS")
print("")
print("FOUND "+dimsizes(brkbad)+" total 'bad' (precip + less than min) times.")
print("")
do f=0,dimsizes(brkbad)-2,1
  dnt = brkbad(f+1)-brkbad(f)
  if dnt.gt.1 then
    nmin = dnt%60
    nhrs = dnt/60
    ;print("")
    ;print("LENGTH OF PERIOD = "+nhrs+" HRS "+nmin+" MIN")
    ;print("PER BEG IDX = "+brkbad(f))
    ;print("PER END IDX = "+brkbad(f+1))
    badbrkprecip(f) = dnt
  end if
end do
npersbrk = dimsizes(ind(.not.ismissing(badbrkprecip)))
print("FOUND "+npersbrk+" GOOD PERIODS.")
print("")
wks4 = gsn_open_wks("png",imgout+"histbrk_"+sprinti("%.3i",toint(mingood)))
bins4 = gsn_histogram(wks4,badbrkprecip,hres)
draw(bins4)
frame(wks4)

; Create a time series plot, showing the number of hours that we'd be averaging from each day of ST output,
; after allowing a minimum period of 4 hours of continuous data.
nptsmon = new((/1440/),float)
nptsbrk = new((/1440/),float)
nptsmon(:) = 0.0
nptsbrk(:) = 0.0
do t=0,1439,1

  ; MONSOON
  tmpm = new((/ndayssoon/),float)
  tmpm@_FillValue = bf
  do nd=0,ndayssoon-1,1
    tmpm(nd) = monflag(((t)+(nd*1440)))
  end do 
  monval = ind(.not.ismissing(tmpm))
  if dimsizes(monval).gt.1 .and. .not.ismissing(monval(0)) then
    nptsmon(t) = (dimsizes(tmpm)-sum(tmpm(monval)))
  else
    delete([/monval,tmpm/])
  end if

  ; BREAK
  tmpb = new((/ndaysbreak/),float)
  tmpb@_FillValue = bf
  do nd=0,ndaysbreak-1,1
    tmpb(nd) = brkflag(((t)+(nd*1440)))
  end do
  brkval = ind(.not.ismissing(tmpb))
  if dimsizes(brkval).gt.1 .and. .not.ismissing(brkval(0)) then
    nptsbrk(t) = (dimsizes(tmpb)-sum(tmpb(brkval)))
  else
    delete([/brkval,tmpb/])
  end if
end do

; Plot two panel plot showing number of times used in ST average for break/monsoon (composite essentially)
wksc = gsn_open_wks("png",imgout+"nst_"+sprinti("%.3i",toint(mingood)))
gsn_define_colormap(wksc,"rainbow+gray")
res = True
res@gsnDraw = False
res@gsnFrame = False
res@vpYF = 0.90
res@vpXF = 0.15
res@vpHeightF = 0.32
res@vpWidthF = 0.75
res@trYMaxF = ndayssoon
res@trYMinF = 0
res@tiYAxisString = "Monsoon Days"
res@tiYAxisFontHeightF = 0.012
res@xyLineThicknessF = 2

xvals = ispan(0,1439,1)
p1 = gsn_csm_y(wksc,nptsmon,res)
draw(p1)

res@vpYF = 0.45
res@trYMaxF = ndaysbreak
res@tiYAxisString = "Break Days"
res@tiYAxisFontHeightF = 0.012
p2 = gsn_csm_y(wksc,nptsbrk,res)
draw(p2)
frame(wksc)

;========================
; USE CODE BELOW TO APPLY THE NEW MASK TO THE DATA
;
; Store the final pass after filtering out "noisy" periods of good data less than threshold
;rainpass2 = new((/dimsizes(raindbz)/),float)
;rainpass2@_FillValue = bf
;omegpass2 = new((/dimsizes(omean)/),float)
;omegpass2@_FillValue = bf
;do nt=0,ntimes-1,1
;  if precipflag(nt).gt.0 then
;    rainpass2(nt,:) = bf
;    omegpass2(nt,:) = bf
;  else
;    rainpass2(nt,:) = raindbz(nt,:)
;    omegpass2(nt,:) = omean(nt,:)
;  end if
;end do
