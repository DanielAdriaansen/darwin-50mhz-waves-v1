;
; File: plot_comp_panel.ncl
;
; Author: D. Adriaansen
;
; Date: 29 Apr 2016
;
; Purpose:
;
; Notes:
;_________________________________________________________________________________

;################################### User Config ##############################;

; Path to input netCDF files
ncpath = "/d1/dadriaan/paper/data/masked"

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
vers = "raw"

; Path to image output
imgout = "/d1/dadriaan/paper/images/comp_"+vers

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
if vers.eq."raw" then
  w = fin[:]->omean
  u = fin[:]->umean
  v = fin[:]->vmean
else
  w = fin[:]->omegpass1
  u = fin[:]->zonalpass1
  v = fin[:]->meridpass1
end if


