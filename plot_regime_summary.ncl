;
; File: plot_regime_summary.ncl
;
; Author: D. Adriaansen
;
; Date: 08 Jun 2016
;
; Purpose: plot two three panel plots of time/height with number of valid data points at each time for
;          all days in the regime for each 50 MHz wind component
;
; Notes: this is intended to show the data after it's been filtered for precip, bad data, and a minimum number of periods
;
;____________________________________________________________________________________________________________________

;##################################### User Config #########################################;

; Path to input netCDF files
ncpath = "/d1/dadriaan/paper/data/c3/maskedminbad"

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

; Variable to make sure we've figured out how to slice the data
foundoffset = False

; What wind component?
windcomp = "w"

; Break or monsoon?
bm = "break"

; What version of the data? (raw, data, prime)?
; raw = raw data
; data = raw data with precip mask applied
; prime = perturbation wind
vers = "prime"

; Path to image output
imgout = "/d1/dadriaan/paper/images"

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

;####################################################################################;

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
;omean = fin[:]->omean
;zmean = fin[:]->zmean
;mmean = fin[:]->mmean
omegpass1 = fin[:]->omegpass2
zonalpass1 = fin[:]->zonalpass2
meridpass1 = fin[:]->meridpass2

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

; Subset begin and end for collecting data
if bm.eq."monsoon" then
  sub_beg = begmonsoon
  sub_end = endmonsoon
else
  sub_beg = begbreak
  sub_end = endbreak
end if
print("CREATING SUMMARY PLOTS FOR DATES FROM:")
print(""+prettyPrint(unix_time(sub_beg),"ymdhms"))
print("TO")
print(""+prettyPrint(unix_time(sub_end),"ymdhms"))
print("")

; Create 2D variable to hold the number of points at each time in the 24-hr period representing each regime
wnpts = new((/1440,numhgts/),float)
wnpts@_FillValue = bf
wnpts!0 = "time"
wnpts!1 = "height"
unpts = new((/1440,numhgts/),float)
unpts@_FillValue = bf
unpts!0 = "time"
unpts!1 = "height"
vnpts = new((/1440,numhgts/),float)
vnpts@_FillValue = bf
vnpts!0 = "time"
vnpts!1 = "height"

; Now for each time, collect all the points at each height that are valid based on the masking that occurred
do t=0,1439,1
  ; At each time, loop over every height and collect data forward in time matching same minutes
  do z=0,numhgts-1,1
    ; Create a temporary array to hold the minutes that we're at from each day in the period (break/monsoon)
    if bm.eq."monsoon" then
      tmpw = new((/ndayssoon/),float)
      tmpu = new((/ndayssoon/),float)
      tmpv = new((/ndayssoon/),float)
      daysinloop = ndayssoon
    else
      tmpw = new((/ndaysbreak/),float)
      tmpu = new((/ndaysbreak/),float)
      tmpv = new((/ndaysbreak/),float)
      daysinloop = ndaysbreak
    end if
    
    ; Loop over every day in the period and collect the relevant data
    do nd=0,daysinloop-1,1
      if bm.eq."monsoon" then
        tmpw(nd) = omegpass1(((begmonsoon+t)+(nd*1440)),z)
        tmpu(nd) = zonalpass1(((begmonsoon+t)+(nd*1440)),z)
        tmpv(nd) = meridpass1(((begmonsoon+t)+(nd*1440)),z)
      else
        tmpw(nd) = omegpass1(((begbreak+t)+(nd*1440)),z)
        tmpu(nd) = zonalpass1(((begbreak+t)+(nd*1440)),z)
        tmpv(nd) = meridpass1(((begbreak+t)+(nd*1440)),z)
      end if
    end do

    ; Determine the amount of good data at this height for each minute in the 24 hour period
    wvalid = ind(.not.ismissing(tmpw))
    uvalid = ind(.not.ismissing(tmpu))
    vvalid = ind(.not.ismissing(tmpv))
 
    ; Store the count
    if dimsizes(wvalid).gt.1 .and. .not.ismissing(wvalid(0)) then
      wnpts(t,z) = dimsizes(wvalid)
      delete([/wvalid,tmpw/])
    else
      delete([/wvalid,tmpw/])
    end if
    if dimsizes(uvalid).gt.1 .and. .not.ismissing(uvalid(0)) then
      unpts(t,z) = dimsizes(uvalid)
      delete([/uvalid,tmpu/])
    else
      delete([/uvalid,tmpu/])
    end if
    if dimsizes(vvalid).gt.1 .and. .not.ismissing(vvalid(0)) then
      vnpts(t,z) = dimsizes(vvalid)
      delete([/vvalid,tmpv/])
    else
      delete([/vvalid,tmpv/])
    end if
  end do
end do

printVarSummary(wnpts)
printVarSummary(unpts)
printVarSummary(vnpts)

printMinMax(wnpts,0)
printMinMax(unpts,0)
printMinMax(vnpts,0)

; Open a workstation for the plot
wks = gsn_open_wks("png",bm+"_regime_summary")

; Resources
res = True
res@gsnMaximize = False
res@gsnFrame = False
res@gsnDraw = False
res@cnFillMode = "RasterFill"
res@cnFillOn = True
res@cnLineLabelsOn = False
res@cnLinesOn = False
res@vpYF = 0.95
res@vpHeightF = 0.25
res@vpWidthF = 0.7
res@vpXF = 0.15
res@cnLevelSelectionMode = "ExplicitLevels"
res@lbLabelBarOn = True
res@cnInfoLabelOn = False
res@tmXTOn = False
res@tmYROn = False
res@gsnLeftStringFontHeightF = 0.012
res@gsnRightStringFontHeightF = 0.012
if bm.eq."monsoon" then
  res@gsnLeftString = prettyPrint(unix_time(begmonsoon),"ymdhm")
  res@gsnRightString = prettyPrint(unix_time(endmonsoon),"ymdhm")
else
  res@gsnLeftString = prettyPrint(unix_time(begbreak),"ymdhm")
  res@gsnRightString = prettyPrint(unix_time(endbreak),"ymdhm")
end if
res@tiXAxisFontHeightF = 0.01
res@tmXBLabelFontHeightF = 0.01
res@tmYLMode = "Explicit"
res@tmYLValues = (/0,20,40,60,80,100,120/)
res@tmYLLabels = tostring(res@tmYLValues/10)
res@tmYLLabelFontHeightF = 0.01
res@tmXBMode = "Explicit"
res@tmXBValues = (/0,120,240,360,480,600,720,840,960,1080,1200,1320,1440/)
res@tmXBLabels = (/"1130","1330","1530","1730","1930","2130","2330","0130","0330","0530","0730","0930","1130"/)
res@tmXBLabelFontHeightF = 0.01
res@trXMaxF = 1440
res@trYMinF = 0.0
res@tiYAxisString = "Height (km)"
res@tiYAxisFontHeightF = 0.01
if bm.eq."monsoon" then
  res@cnLevels = ispan(0,ndayssoon-1,1)
  labs = new((/dimsizes(res@cnLevels)+2/),string)
  labs(0) = ""
  labs(1:ndayssoon) = tostring(res@cnLevels)
  labs(dimsizes(labs)-1) = ""
  res@lbLabelStrings = labs
  ;tres@cnLevels = (/0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18/)
  ;tres@lbLabelStrings = (/"","0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18",""/)
else
  res@cnLevels = ispan(0,ndaysbreak-1,1)
  labs = new((/dimsizes(res@cnLevels)+2/),string)
  labs(0) = ""
  labs(1:ndaysbreak) = tostring(res@cnLevels)
  labs(dimsizes(labs)-1) = ""
  res@lbLabelStrings = labs
  ;tres@cnLevels = (/0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21/)
  ;tres@lbLabelStrings = (/"","0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21",""/)
end if
res@cnFillPalette = "ncl_default"
res@lbOrientation = "vertical"
res@lbLeftMarginF = 0.30
res@lbRightMarginF = 0.0
res@lbTitleFontHeightF = 0.010
res@lbLabelFontHeightF = 0.008
res@lbTitlePosition = "Right"
res@lbTitleDirection = "Across"
res@lbTitleAngleF = 90.0
res@lbTitleOffsetF = 0.3
res@lbTitleString = "W Count"
res@pmLabelBarOrthogonalPosF = 0.0
res@pmLabelBarWidthF = 0.07
res@pmLabelBarSide = "Right"
res@lbLabelAlignment = "ExternalEdges"
res@tiMainOffsetYF = -0.016

p1 = gsn_csm_contour(wks,wnpts(height|0:120,time|:),res)
draw(p1)
res@vpYF = 0.64
res@gsnLeftString = ""
res@gsnRightString = ""
res@tiMainString = ""
res@lbTitleString = "U Count"
p2 = gsn_csm_contour(wks,unpts(height|0:120,time|:),res)
draw(p2)
res@vpYF = 0.33
res@tiXAxisString = "Local Time (Darwin)"
res@lbTitleString = "V Count"
p3 = gsn_csm_contour(wks,vnpts(height|0:120,time|:),res)
draw(p3)
frame(wks)
