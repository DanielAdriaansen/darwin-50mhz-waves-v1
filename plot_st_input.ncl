;
; File: plot_data_panel.ncl
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
;ncpath = "/d1/dadriaan/paper/data/maskedmin"
ncpath = "/d1/dadriaan/paper/data/maskedminbad"

; Path to image output
;imgout = "/d1/dadriaan/paper/images/stdata"
imgout = "/d1/dadriaan/paper/images/stdata2"

; What hour is the beginning of a day? In Darwin, we will use 02Z to 02Z, or 1130-1130 local time.
beghr = 2

; What day number do we want to look at? Use any negative number to indicate we want all days
if .not.isvar("daynum") then
  daynum = 23
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
omegpass1 = fin[:]->omegpass2
zonalpass1 = fin[:]->zonalpass2
meridpass1 = fin[:]->meridpass2
rainpass1 = fin[:]->rainpass2

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
  wks = gsn_open_wks("png",imgout+"/stdata_"+sprinti("%.2i",beghr)+"Z")
else
  wks = gsn_open_wks("png",imgout+"/stdata_"+sprinti("%.2i",beghr)+"Z_"+sprinti("%.2i",daynum))
end if

; Resources
pres = True
pres@gsnMaximize = False
pres@gsnFrame = False
pres@gsnDraw = False
pres@cnFillMode = "RasterFill"
pres@cnFillOn = True
pres@cnLineLabelsOn = False
pres@cnLinesOn = False
;pres@vpYF = 0.95 ; THREE TOP
;pres@vpHeightF = 0.25 ; THREE PANELS
pres@vpYF = 0.97 ; FOUR TOP
pres@vpHeightF = 0.18 ; FOUR PANELS
pres@vpWidthF = 0.7
pres@vpXF = 0.15
pres@cnLevelSelectionMode = "ExplicitLevels"
pres@lbLabelBarOn = True
pres@cnInfoLabelOn = False
pres@tmXTOn = False
pres@tmYROn = False
pres@gsnLeftStringFontHeightF = 0.012
pres@gsnRightStringFontHeightF = 0.012
pres@gsnLeftString = sstring
pres@gsnRightString = estring
pres@tiXAxisFontHeightF = 0.01
if daynum.lt.0 then
  pres@tmXBMode = "Explicit"
  pres@tmXBValues = daybeg
  pres@tmXBLabels = tostring(doy)
  pres@tmXBLabelFontHeightF = 0.01
  pres@tmXBLabelAngleF = -90.0
  pres@tmXBLabelJust = "CenterCenter"
else
  pres@tmXBLabelFontHeightF = 0.01
end if
pres@tmYLMode = "Explicit"
pres@tmYLValues = (/0,20,40,60,80,100,120/)
pres@tmYLLabels = tostring(pres@tmYLValues/10)
pres@tmYLLabelFontHeightF = 0.01
pres@tmXBMode = "Explicit"
pres@tmXBValues = (/0,120,240,360,480,600,720,840,960,1080,1200,1320,1440/)
pres@tmXBLabels = (/"1130","1330","1530","1730","1930","2130","2330","0130","0330","0530","0730","0930","1130"/)
pres@tmXBLabelFontHeightF = 0.01
pres@trYMinF = 0.0
pres@tiYAxisString = "Height (km)"
pres@tiYAxisFontHeightF = 0.01
;pres@cnLevels = fspan(-0.5,0.5,20)
pres@cnLevels = (/-0.5,-0.4,-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3,0.4,0.5/)
pres@cnFillPalette = "ncl_default"
pres@lbLabelStrings = (/"","-0.5","-0.4","-0.3","-0.2","-0.1","0.0","0.1","0.2","0.3","0.4","0.5",""/)
pres@lbOrientation = "vertical"
pres@lbLeftMarginF = 0.30
pres@lbRightMarginF = 0.0
pres@lbTitleFontHeightF = 0.010
pres@lbLabelFontHeightF = 0.008
pres@lbTitlePosition = "Right"
pres@lbTitleDirection = "Across"
pres@lbTitleAngleF = 90.0
pres@lbTitleOffsetF = 0.3
pres@pmLabelBarOrthogonalPosF = 0.0
pres@pmLabelBarWidthF = 0.07
pres@pmLabelBarSide = "Right"
pres@lbLabelAlignment = "BoxCenters"
pres@lbTitleString = "Type"
pres@lbLabelAlignment = "ExternalEdges"
pres@lbTitleString = "50 MHz w (m/s)" 

; PLOTS
;p12 = gsn_csm_contour(wks,omean(height|0:120,time|beg_ind:end_ind),pres)
p12 = gsn_csm_contour(wks,omegpass1(height|0:120,time|beg_ind:end_ind),pres)
draw(p12)
;pres@vpYF = 0.64 ; THREE MIDDLE
pres@vpYF = 0.73 ; FOUR MIDDLE TOP
pres@gsnLeftString = ""
pres@gsnRightString = ""
delete(pres@cnLevels)
delete(pres@lbLabelStrings)
;pres@cnLevels = fspan(-20,20,20)
pres@cnLevels = (/-10,-8,-6,-4,-2,0,2,4,6,8,10/)
pres@lbLabelStrings = (/"","-10","-8","-6","-4","-2","0","2","4","6","8","10",""/)
pres@lbTitleString = "50 MHz u (m/s)"
;p22 = gsn_csm_contour(wks,zmean(height|0:120,time|beg_ind:end_ind),pres)
p22 = gsn_csm_contour(wks,zonalpass1(height|0:120,time|beg_ind:end_ind),pres)
draw(p22)
;pres@vpYF = 0.33 ; THREE BOTTOM
pres@vpYF = 0.49 ; FOUR MIDDLE BOTTOM
pres@lbTitleString = "50 MHz v (m/s)"
;p32 = gsn_csm_contour(wks,mmean(height|0:120,time|beg_ind:end_ind),pres)
p32 = gsn_csm_contour(wks,meridpass1(height|0:120,time|beg_ind:end_ind),pres)
draw(p32)
pres@vpYF = 0.25 ; FOUR BOTTOM
;pres@tiXAxisString = "Time (minutes)"
if daynum.lt.0 then
  pres@tiXAxisString = "Day of year"
else
  pres@tiXAxisString = "Local Time (Darwin)"
end if
delete(pres@cnLevels)
delete(pres@lbLabelStrings)
pres@cnLevels = (/0.0,5.0,10.0,15.0,20.0,25.0,30.0,35.0,40.0,45.0,50.0/)
pres@lbLabelStrings = (/"","0","5","10","15","20","25","30","35","40","45","50",""/)
pres@lbTitleString = "920 MHz (dBZ)"
;p42 = gsn_csm_contour(wks,raindbz(height|0:120,time|beg_ind:end_ind),pres)
p42 = gsn_csm_contour(wks,rainpass1(height|0:120,time|beg_ind:end_ind),pres)
draw(p42)
frame(wks)
