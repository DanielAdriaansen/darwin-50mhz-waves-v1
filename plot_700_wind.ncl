;
; File: plot_700_wind.ncl
;
; Author: D. Adriaansen
;
; Date: 05 May 2016
;
; Purpose: Read in 50 MHz wind data, find 700 hPa (approx), take a 4 day running mean and plot the time series
;
; Notes: Monsoon vs. break methodology taken from TWP-ICE BAMS paper
;_________________________________________________________________________________

;################################### User Config ##################################;

; Path to input netCDF files
;ncpath = "/d1/dadriaan/paper/data/masked"
ncpath = "/d1/dadriaan/paper/data/c2/masked"

; Path to image output
imgout = "/d1/dadriaan/paper/images/raw"

; What hour is the beginning of a day? In Darwin, we will use 02Z to 02Z, or 1130-1130 local time.
beghr = 2

; What day number do we want to look at? Use any negative number to inicate we want all days
if .not.isvar("daynum") then
  daynum = -1
end if

; Variable to make sure we've figure out how to slice the data
foundoffset = False

; What level for monsoon definition? Level 29 is approx 3000m, or 700 hPa. Note, level "0" = 100m.
mlev = 29

; Print configuration for user:
print("")
print("###############################")
print("")
print("DAY BEGINS AT: "+sprinti("%02i",tointeger(beghr))+"Z")
if daynum.lt.0 then
  print("PLOTTING DATA FOR ALL DAYS")
else
  print("PLOTTING DATA FOR DAY NUMBER: "+daynum)
end if
print("SELECTING DATA FROM LEVEL "+mlev)
print("")
print("################################")
print("")

;####################################################################################;

; NCL Libraries and functions
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
omean = fin[:]->omean
zmean = fin[:]->zmean
mmean = fin[:]->mmean

; Print info for user
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

; Before we plot the mean wind, let's take a four day running mean
windseries = zmean(:,mlev)
windseries@_FillValue = bf
smoothwind = new((/dimsizes(windseries)/),float)
offset = round(((1440*4)/2),3)
print("RUNNING MEAN OFFSET = "+offset)
print("")
do i=0,dimsizes(smoothwind)-1,1
  ;print("i = "+i)
  ; If we don't have a full 1440*4 number of points to use, just retain the original data
  if i.lt.offset then
    ;print("TOO CLOSE TO BEG. RETAINING DATA")
    smoothwind(i) = windseries(i)
    continue
  else
    if i.gt.(dimsizes(smoothwind)-1)-offset then
      ;print("TOO CLOSE TO END. SETTING TO MISSING")
      smoothwind(i) = bf
    else
      ; Otherwise do the 4 day mean (+/- from current point)
      ;print("ATTEMPTING MEAN")
      subset = windseries(i-offset:i+offset)
      subset@_FillValue = bf
      nomiss = ind(.not.ismissing(subset))
      ;print("NUM GOOD ="+dimsizes(nomiss))
      if dimsizes(nomiss).gt.0 .and. .not.ismissing(nomiss(0)) then
        ;print("FOUND SOME GOOD")
        avgwind = subset(nomiss)
      else
        ;print("ALL MISSING")
        smoothwind(i) = bf
        delete(subset)
        delete(nomiss)
        continue
      end if
      avgwind = subset(nomiss)
      if dimsizes(avgwind).lt.offset then
        smoothwind(i) = bf
        ;print("NOT ENOUGH POINTS")
        delete(subset)
        delete(avgwind)
        delete(nomiss)
        continue
      else
        ;print("USING "+dimsizes(avgwind)+" POINTS FOR MEAN")
        smoothwind(i) = avg(avgwind)
      end if
      delete(subset)
      delete(nomiss)
      delete(avgwind)
    end if
  end if
end do

;----------------------------------------- PLOTTING FOR PROFILER MEAN WIND
zres = True
if daynum.lt.0 then
  zwks = gsn_open_wks("png","zpanel_"+sprinti("%.2i",beghr)+"Z")
else
  zwks = gsn_open_wks("png","zpanel_"+sprinti("%.2i",beghr)+"Z_"+sprinti("%.2i",daynum))
end if
gsn_define_colormap(zwks,"precip4_11lev")
zres@gsnMaximize = False
zres@gsnDraw = False
zres@gsnFrame = False
;zres@trXMaxF = dimsizes(zmean(beg_ind:end_ind,mlev))-1
zres@trXMaxF = dimsizes(smoothwind(beg_ind:end_ind))-1
zres@trXMinF = 0
zres@trYMaxF = 40.0
zres@trYMinF = -40.0
zres@tiYAxisString = "Zonal Wind (m/s) @ 3000m (approx. 700 hPa)"
if daynum.lt.0 then
  zres@tiXAxisString = "Day of year"
else
  zres@tiXAxisString = "Time (minutes)"
end if
zres@tiXAxisFontHeightF = 0.012
zres@tiYAxisFontHeightF = 0.012
zres@tiMainFontHeightF = 0.014
zres@tiMainString = "Positive values are eastward"
zres@tmXBLabelFontHeightF = 0.012
zres@tmYLLabelFontHeightF = 0.012
zres@gsnLeftStringFontHeightF = 0.012
zres@gsnRightStringFontHeightF = 0.012
zres@gsnLeftString = sstring
zres@gsnRightString = estring
zres@xyLineThicknessF = 4
zres@xyLineColor = "black"
if daynum.lt.0 then
  zres@tmXBMode = "Explicit"
  zres@tmXBValues = daybeg
  zres@tmXBLabels = tostring(doy)
  zres@tmXBLabelFontHeightF = 0.01
  zres@tmXBLabelAngleF = -90.0
  zres@tmXBLabelJust = "CenterCenter"
end if
zres@vpWidthF = 0.7

zres@tmXTOn = False
zres@tmYROn = False

zres@tmXMajorGrid = True
zres@tmXMajorGridLineColor = "black"
zres@tmXMajorGridLineDashPattern = 0
zres@tmXMajorGridThicknessF = 1

zres@tmYMajorGrid = True
zres@tmYMajorGridLineColor = "black"
zres@tmYMajorGridLineDashPattern = 0
zres@tmYMajorGridThicknessF = 1

zres@vpYF = 0.90
zres@vpXF = 0.15
zres@vpHeightF = 0.32
zres@vpWidthF = 0.75

zres@tfPolyDrawOrder = "PreDraw"

; Resources for polygons/lines
pgres = True
pgres@gsFillColor = 3

;pgres@gsLineColor = "red"
;pgres@gsLineThicknessF = 2
;liney = new((dimsizes(ispan(toint(zres@trYMinF),toint(zres@trYMaxF),2))),float)
;liney(:) = ispan(toint(zres@trYMinF),toint(zres@trYMaxF),2)
;linex = new((dimsizes(liney)),float)
;linex(:) = daybeg(37)

; Break, days 37 - end
;pgbx = (/daybeg(36),daybeg(36),daybeg(dimsizes(daybeg)-1),daybeg(dimsizes(daybeg)-1),daybeg(36)/)
pgbx = (/daybeg(36),daybeg(36),zres@trXMaxF,zres@trXMaxF,daybeg(36)/)
pgby = (/zres@trYMinF,zres@trYMaxF,zres@trYMaxF,zres@trYMinF,zres@trYMinF/)

; Monsoon, days 16 - 33
;pgmx = (/daybeg(15),daybeg(15),daybeg(32),daybeg(32),daybeg(15)/)
pgmx = (/daybeg(13),daybeg(13),daybeg(33),daybeg(33),daybeg(13)/)
pgmy = (/zres@trYMinF,zres@trYMaxF,zres@trYMaxF,zres@trYMinF,zres@trYMinF/)

p1 = gsn_csm_y(zwks,windseries(beg_ind:end_ind),zres)
dumb1 = gsn_add_polygon(zwks,p1,pgbx,pgby,pgres)
pgres@gsFillColor = 8
dumm1 = gsn_add_polygon(zwks,p1,pgmx,pgmy,pgres)
draw(p1)

zres@vpYF = 0.45
zres@tiMainString = "4 Day Running Mean"
pgres@gsFillColor = 3
p2 = gsn_csm_y(zwks,smoothwind(beg_ind:end_ind),zres)
dumb2 = gsn_add_polygon(zwks,p2,pgbx,pgby,pgres)
pgres@gsFillColor = 8
dumm2 = gsn_add_polygon(zwks,p2,pgmx,pgmy,pgres)
;draw(p2)

; Code for filled square marker
mstring = "y"
fontnum = 35
xoffset = 0.0
yoffset = 0.0
ratio   = 1.0
size    = 1.0
angle   = 0.0
new_index = NhlNewMarker(zwks, mstring, fontnum, xoffset, yoffset, ratio, size, angle)

lgres = True
lgres@lgAutoManage = False
lgres@lgMarkerColors = 8
lgres@lgItemTypes = "Markers"
lgres@lgMarkerIndexes = new_index
lgres@vpWidthF = 0.15
lgres@lgPerimOn = False
lgres@lgMonoMarkerIndex = False
lgres@lgMonoItemType = False
lgres@lgMonoMarkerColor = False
lgres@lgLabelOffsetF = -0.12
lgres@lgLabelFontHeightF = 0.08
lgres@vpHeightF = 0.1
lbid1 = gsn_create_legend(zwks,1,(/" Monsoon"/),lgres)
delete(lgres)

lgres = True
lgres@lgAutoManage = False
lgres@lgMarkerColors = 3
lgres@lgItemTypes = "Markers"
lgres@lgMarkerIndexes = new_index
lgres@vpWidthF = 0.15
lgres@lgPerimOn = False
lgres@lgMonoMarkerIndex = False
lgres@lgMonoItemType = False
lgres@lgMonoMarkerColor = False
lgres@lgLabelOffsetF = -0.15
lgres@lgLabelFontHeightF = 0.06
lgres@vpHeightF = 0.1
lbid2 = gsn_create_legend(zwks,1,(/" Break"/),lgres)

amres = True
amres@amParallelPosF = -0.15
amres@amOrthogonalPosF = 0.78
annoid1 = gsn_add_annotation(p2,lbid1,amres)
amres@amParallelPosF = 0.15
annoid2 = gsn_add_annotation(p2,lbid2,amres)
draw(p2)
;drawNDCGrid(zwks)
frame(zwks)

print("PROCESSING COMPLETE.")
print("")
