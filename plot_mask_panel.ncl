;
; File: plot_mask_panel.ncl
;
; Author: D. Adriaansen
;
; Date: 29 Apr 2016
;
; Purpose: Read in the profiler data and create three panel plots of the data.
;          PANEL 1: 50 MHz w wind (raw)
;          PANEL 2: 50 MHz u wind (raw)
;          PANEL 3: 50 MHz v wind (raw)
;          PANEL 4: 920 MHz DBz (raw)
;
; Notes:
;__________________________________________________________________________________

;############################## User Config ######################################;

; Path to input netCDF files
ncpath = "/d1/dadriaan/paper/data/c3/masked"

; Path to image output
imgout = "/d1/dadriaan/paper/images/mask/c3"

; What hour is the beginning of a day? In Darwin, we will use 02Z to 02Z, or 1130-1130 local time.
beghr = 2

; What day number do we want to look at? Use any negative number to indicate we want all days
if .not.isvar("daynum") then
  daynum = -1
end if

; Variable to make sure we've figured out how to slice the data
foundoffset = False

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

;##################################################################################;

; NCL Libraries
load "/home/dadriaan/cvs/pat/dadriaan/scripts/ncl/load_ncl_libs.ncl"
load "/home/dadriaan/cvs/pat/dadriaan/scripts/ncl/ncl_usr_funcs.ncl"

; Open the files
allfiles = systemfunc("ls "+ncpath+"/*2006*.nc")
fin = addfiles(allfiles,"r")
ListSetType(fin,"cat")

; Read the data
bf = fin[0]->bf
pagl = fin[0]->pagl
unix_time = fin[:]->unix_time
mask_w = fin[:]->mask_w
mask_u = fin[:]->mask_u
mask_v = fin[:]->mask_v
mask_dbz = fin[:]->mask_dbz

; Info for user
nhrsdata = dimsizes(unix_time)/60
ndaysdata = dimsizes(unix_time)/1440
print("PROCESSING "+ndaysdata+" DAYS OF DATA ("+nhrsdata+" HOURS)")
print("TOTAL MINUTES = "+dimsizes(unix_time))
print("")

; Determine the number of times and heights
ntimes = dimsizes(unix_time)
numhgts = dimsizes(pagl)

; Compute the correct start and end index for the plot
; CASE: All days with offset start time (e.g. 02Z)
if daynum.lt.0 .and. beghr.gt.0 then
  beg_ind = (beghr*60)
  end_ind = ntimes-1
  print("CASE 1")
  print("beg_ind = "+beg_ind)
  print("end_ind = "+end_ind)
  print("")
  foundoffset = True
end if
; CASE: All days starting at the beginning of the data (no offset start time)
if daynum.lt.0 .and. beghr.lt.1 then
  beg_ind = 0
  end_ind = ntimes-1
  print("CASE 2")
  print("beg_ind = "+beg_ind)
  print("end_ind = "+end_ind)
  print("")
  foundoffset = True
end if
; CASE: Specific day with offset start time (e.g. 02Z)
if daynum.ge.0 .and. beghr.gt.0 then
  beg_ind = (beghr*60)+(1440*daynum)
  end_ind = beg_ind+1440
  print("CASE 3")
  print("beg_ind = "+beg_ind)
  print("end_ind = "+end_ind)
  print("")
  foundoffset = True
end if
if .not.foundoffset then
  print("FATAL! UNSURE WHAT TO DO WITH daynum AND beghr COMBO.")
  print("REVIEW CONFIGURATION AND RUN AGAIN.")
  print("")
  exit()
end if

; Figure out ymd, hms for plot titles
sstring = prettyPrint(unix_time(beg_ind),"ymdhms")
estring = prettyPrint(unix_time(end_ind),"ymdhms")
print(""+sstring)
print(""+estring)
print("")

; Convert all the UNIX times to calendar objects
allunix = cd_calendar(unix_time,0)

; Create a vector of values starting at the beginning index every 1 minute until the end
xvals = ispan(beghr*60,ntimes-1,1)

; Find the beginning of each day, using the beghr offset
daybeg = ind((xvals-(beghr*60))%1440.eq.0.0)

; Convert the daybeg indices to a day of year, using the unix_time
doy = day_of_year((/toint(allunix(daybeg,0))/),(/toint(allunix(daybeg,1))/),(/toint(allunix(daybeg,2))/))

; Plot the three panel data plot
if daynum.lt.0 then
  wks = gsn_open_wks("png",imgout+"/mask_"+sprinti("%.2i",beghr)+"Z")
else
  wks = gsn_open_wks("png",imgout+"/mask_"+sprinti("%.2i",beghr)+"Z_"+sprinti("%.2i",daynum))
end if
gsn_define_colormap(wks,"rainbow+gray")

; Resources
mres = True
mres@gsnMaximize = False
mres@gsnFrame = False
mres@gsnDraw = False
mres@cnFillMode = "RasterFill"
mres@cnFillOn = True
mres@cnLineLabelsOn = False
mres@cnLinesOn = False
;mres@vpYF = 0.95 ; THREE TOP
;mres@vpHeightF = 0.25 ; THREE PANELS
mres@vpYF = 0.97 ; FOUR TOP
mres@vpHeightF = 0.18 ; FOUR PANELS
mres@vpWidthF = 0.7
mres@vpXF = 0.15
mres@cnLevelSelectionMode = "ExplicitLevels"
mres@cnLevels = (/1.5,2.5/)
mres@cnFillColors = (/130,180,238/)
mres@lbLabelBarOn = True
mres@cnInfoLabelOn = False
mres@tmXTOn = False
mres@tmYROn = False
mres@gsnLeftStringFontHeightF = 0.012
mres@gsnRightStringFontHeightF = 0.012
mres@gsnLeftString = sstring
mres@gsnRightString = estring
mres@tiXAxisFontHeightF = 0.01
if daynum.lt.0 then
  mres@tmXBMode = "Explicit"
  mres@tmXBValues = daybeg
  mres@tmXBLabels = tostring(doy)
  mres@tmXBLabelFontHeightF = 0.01
  mres@tmXBLabelAngleF = -90.0
  mres@tmXBLabelJust = "CenterCenter"
else
  mres@tmXBLabelFontHeightF = 0.01
end if
mres@tmYLMode = "Explicit"
mres@tmYLValues = (/0,20,40,60,80,100,120/)
mres@tmYLLabels = tostring(mres@tmYLValues/10)
mres@tmYLLabelFontHeightF = 0.01
mres@tmXBMode = "Explicit"
mres@tmXBValues = (/0,120,240,360,480,600,720,840,960,1080,1200,1320,1440/)
mres@tmXBLabels = (/"1130","1330","1530","1730","1930","2130","2330","0130","0330","0530","0730","0930","1130"/)
mres@tmXBLabelFontHeightF = 0.01
mres@trYMinF = 0.0
mres@tiYAxisString = "Height (km)"
mres@tiYAxisFontHeightF = 0.01
mres@lbOrientation = "vertical"
mres@lbLeftMarginF = 0.30
mres@lbRightMarginF = 0.0
mres@lbTitleFontHeightF = 0.010
mres@lbLabelFontHeightF = 0.008
mres@lbTitlePosition = "Right"
mres@lbTitleDirection = "Across"
mres@lbTitleAngleF = 90.0
mres@lbTitleOffsetF = 0.3
mres@pmLabelBarOrthogonalPosF = 0.0
mres@pmLabelBarWidthF = 0.07
mres@pmLabelBarSide = "Right"
mres@lbLabelStrings = (/"GOOD","BAD","PRECIP"/)
mres@lbLabelAngleF = 90.0
mres@lbLabelAlignment = "BoxCenters"
mres@lbTitleString = ""

; PLOTS
p1 = gsn_csm_contour(wks,mask_w(height|0:120,time|beg_ind:end_ind),mres)
draw(p1)
;mres@vpYF = 0.64 ; THREE MIDDLE
mres@vpYF = 0.73 ; FOUR MIDDLE TOP
mres@gsnLeftString = ""
mres@gsnRightString = ""
p2 = gsn_csm_contour(wks,mask_u(height|0:120,time|beg_ind:end_ind),mres)
draw(p2)
;mres@vpYF = 0.33 ; THREE BOTTOM
mres@vpYF = 0.49 ; FOUR MIDDLE BOTTOM
p3 = gsn_csm_contour(wks,mask_v(height|0:120,time|beg_ind:end_ind),mres)
draw(p3)
mres@vpYF = 0.25 ; FOUR BOTTOM
;mres@tiXAxisString = "Time (minutes)"
if daynum.lt.0 then
  mres@tiXAxisString = "Day of year"
else
  mres@tiXAxisString = "Local Time (Darwin)"
end if
p4 = gsn_csm_contour(wks,mask_dbz(height|0:120,time|beg_ind:end_ind),mres)
draw(p4)
frame(wks)
