;
; File: plot_prime_panel.ncl
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
ncpath = "/d1/dadriaan/paper/data/masked"

; Path to image output
imgout = "/d1/dadriaan/paper/images/prime"

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
omegpass1 = fin[:]->omegpass1
zonalpass1 = fin[:]->zonalpass1
meridpass1 = fin[:]->meridpass1
rainpass1 = fin[:]->rainpass1

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

; Take the mean at each height for each wind component
meanw = new((/numhgts/),float)
sdevw = new((/numhgts/),float)
nptsw = new((/numhgts/),float)
meanu = new((/numhgts/),float)
sdevu = new((/numhgts/),float)
nptsu = new((/numhgts/),float)
meanv = new((/numhgts/),float)
sdevv = new((/numhgts/),float)
nptsv = new((/numhgts/),float)
do z=0,numhgts-1,1
  ;print("CALCULATING MEANS FOR Z = "+pagl(z))
  subw = omegpass1(beg_ind:end_ind,z)
  subu = zonalpass1(beg_ind:end_ind,z)
  subv = meridpass1(beg_ind:end_ind,z)
  subu@_FillValue = bf
  subw@_FillValue = bf
  subv@_FillValue = bf
  goodw = ind(.not.ismissing(subw))
  goodu = ind(.not.ismissing(subu))
  goodv = ind(.not.ismissing(subv))
  if ismissing(goodw(0)) then
    meanw(z) = bf
    sdevw(z) = bf
    nptsw(z) = bf
  else
    ;print("CALCULATING MEAN WWIND WITH "+dimsizes(goodw)+" POINTS.")
    meanw(z) = avg(subw(goodw))
    sdevw(z) = stddev(subw(goodw))
    nptsw(z) = dimsizes(goodw)
  end if
  ;print("meanw = "+meanw(z))
  ;print("sdevw = "+sdevw(z))
  if ismissing(goodu(0)) then
    meanu(z) = bf
    sdevu(z) = bf
    nptsu(z) = bf
  else
    ;print("CALCULATING MEAN UWIND WITH "+dimsizes(goodu)+" POINTS.")
    meanu(z) = avg(subu(goodu))
    sdevu(z) = stddev(subu(goodu))
    nptsu(z) = dimsizes(goodu)
  end if
  ;print("meanu = "+meanu(z))
  ;print("sdevu = "+sdevu(z))
  if ismissing(goodv(0)) then
    meanv(z) = bf
    sdevv(z) = bf
    nptsv(z) = bf
  else
    ;print("CALCULATING MEAN VWIND WITH "+dimsizes(goodv)+" POINTS.")
    meanv(z) = avg(subv(goodv))
    sdevv(z) = stddev(subv(goodv))
    nptsv(z) = dimsizes(goodv)
  end if
  ;print("meanv = "+meanv(z))
  ;print("sdevv = "+sdevv(z))
  ;print("")

  ; Clean up
  delete([/goodw,goodv,goodu/])
end do

; Define new matrices to hold the means
wmean = new((/dimsizes(omegpass1)/),float)
umean = new((/dimsizes(zonalpass1)/),float)
vmean = new((/dimsizes(meridpass1)/),float)

; Fill the means
do nz=0,numhgts-1,1
  wmean(:,nz) = meanw(nz)
  umean(:,nz) = meanu(nz)
  vmean(:,nz) = meanv(nz)
end do

; Subtract off the mean
wturb = omegpass1-wmean
uturb = zonalpass1-umean
vturb = meridpass1-vmean

; Name dimensions
wturb!0 = "time"
wturb!1 = "height"
vturb!0 = "time"
vturb!1 = "height"
uturb!0 = "time"
uturb!1 = "height"

print("PERTURBATION COMPLETE.")
print("")

; Plot the three panel data plot
if daynum.lt.0 then
  wks = gsn_open_wks("png",imgout+"/prime_"+sprinti("%.2i",beghr)+"Z")
else
  wks = gsn_open_wks("png",imgout+"/prime_"+sprinti("%.2i",beghr)+"Z_"+sprinti("%.2i",daynum))
end if

; Resources
;gsn_define_colormap(pwks,"rainbow+gray")
tres = True
tres@gsnMaximize = False
tres@gsnFrame = False
tres@gsnDraw = False
tres@cnFillMode = "RasterFill"
tres@cnFillOn = True
tres@cnLineLabelsOn = False
tres@cnLinesOn = False
;tres@vpYF = 0.95 ; THREE TOP
;tres@vpHeightF = 0.25 ; THREE PANELS
tres@vpYF = 0.97 ; FOUR TOP
tres@vpHeightF = 0.18 ; FOUR PANELS
tres@vpWidthF = 0.7
tres@vpXF = 0.15
tres@cnLevelSelectionMode = "ExplicitLevels"
;pres@cnLevelSelectionMode = "Automatic"
tres@lbLabelBarOn = True
tres@cnInfoLabelOn = False
tres@tmXTOn = False
tres@tmYROn = False
tres@gsnLeftStringFontHeightF = 0.012
tres@gsnRightStringFontHeightF = 0.012
tres@gsnLeftString = sstring
tres@gsnRightString = estring
tres@tiXAxisFontHeightF = 0.01
if daynum.lt.0 then
  tres@tmXBMode = "Explicit"
  tres@tmXBValues = daybeg
  tres@tmXBLabels = tostring(doy)
  tres@tmXBLabelFontHeightF = 0.01
  tres@tmXBLabelAngleF = -90.0
  tres@tmXBLabelJust = "CenterCenter"
else
  tres@tmXBLabelFontHeightF = 0.01
end if
tres@tmYLMode = "Explicit"
tres@tmYLValues = (/0,20,40,60,80,100,120/)
tres@tmYLLabels = tostring(tres@tmYLValues/10)
tres@tmYLLabelFontHeightF = 0.01
tres@tmXBMode = "Explicit"
tres@tmXBValues = (/0,120,240,360,480,600,720,840,960,1080,1200,1320,1440/)
tres@tmXBLabels = (/"1130","1330","1530","1730","1930","2130","2330","0130","0330","0530","0730","0930","1130"/)
tres@tmXBLabelFontHeightF = 0.01
tres@trYMinF = 0.0
tres@tiYAxisString = "Height (km)"
tres@tiYAxisFontHeightF = 0.01
;tres@cnLevels = fspan(-0.5,0.5,20)
;tres@cnLevels = (/-1.0,-0.8,-0.6,-0.4,-0.2,0.0,0.2,0.4,0.6,0.8,1.0/)
;tres@cnFillPalette = "ncl_default"
;tres@lbLabelStrings = (/"","-1.0","-0.8","-0.6","-0.4","-0.2","0.0","0.2","0.4","0.6","0.8","1.0",""/)
tres@cnLevels = (/-0.5,-0.4,-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3,0.4,0.5/)
tres@cnFillPalette = "ncl_default"
tres@lbLabelStrings = (/"","-0.5","-0.4","-0.3","-0.2","-0.1","0.0","0.1","0.2","0.3","0.4","0.5",""/)
tres@lbOrientation = "vertical"
tres@lbLeftMarginF = 0.30
tres@lbRightMarginF = 0.0
tres@lbTitleFontHeightF = 0.010
tres@lbLabelFontHeightF = 0.008
tres@lbTitlePosition = "Right"
tres@lbTitleDirection = "Across"
tres@lbTitleAngleF = 90.0
tres@lbTitleOffsetF = 0.3
tres@pmLabelBarOrthogonalPosF = 0.0
tres@pmLabelBarWidthF = 0.07
tres@pmLabelBarSide = "Right"
tres@lbLabelAlignment = "ExternalEdges"
tres@lbTitleString = "50 MHz w' (m/s)"

; PLOTS
;t12 = gsn_csm_contour(wks,omean(height|0:120,time|beg_ind:end_ind),tres)
t12 = gsn_csm_contour(wks,wturb(height|0:120,time|beg_ind:end_ind),tres)
draw(t12)
;tres@vpYF = 0.64 ; THREE MIDDLE
tres@vpYF = 0.73 ; FOUR MIDDLE TOP
tres@gsnLeftString = ""
tres@gsnRightString = ""
delete(tres@cnLevels)
delete(tres@lbLabelStrings)
;tres@cnLevels = fspan(-20,20,20)
tres@cnLevels = (/-10,-8,-6,-4,-2,0,2,4,6,8,10/)
tres@lbLabelStrings = (/"","-10","-8","-6","-4","-2","0","2","4","6","8","10",""/)
tres@lbTitleString = "50 MHz u' (m/s)"
;t22 = gsn_csm_contour(wks,zmean(height|0:120,time|beg_ind:end_ind),tres)
t22 = gsn_csm_contour(wks,uturb(height|0:120,time|beg_ind:end_ind),tres)
draw(t22)
;tres@vpYF = 0.33 ; THREE BOTTOM
tres@vpYF = 0.49 ; FOUR MIDDLE BOTTOM
tres@lbTitleString = "50 MHz v' (m/s)"
;t32 = gsn_csm_contour(wks,mmean(height|0:120,time|beg_ind:end_ind),tres)
t32 = gsn_csm_contour(wks,vturb(height|0:120,time|beg_ind:end_ind),tres)
draw(t32)
tres@vpYF = 0.25 ; FOUR BOTTOM
;tres@tiXAxisString = "Time (minutes)"
if daynum.lt.0 then
  tres@tiXAxisString = "Day of year"
else
  tres@tiXAxisString = "Local Time (Darwin)"
end if
delete(tres@cnLevels)
delete(tres@lbLabelStrings)
tres@cnLevels = (/0.0,5.0,10.0,15.0,20.0,25.0,30.0,35.0,40.0,45.0,50.0/)
tres@lbLabelStrings = (/"","0","5","10","15","20","25","30","35","40","45","50",""/)
tres@lbTitleString = "920 MHz (dBZ)"
;t42 = gsn_csm_contour(wks,raindbz(height|0:120,time|beg_ind:end_ind),tres)
t42 = gsn_csm_contour(wks,rainpass1(height|0:120,time|beg_ind:end_ind),tres)
draw(t42)
frame(wks)
