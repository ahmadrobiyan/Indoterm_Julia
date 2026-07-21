echo on
REM this BAT runs TABLO and LTG for %1.TAB only IF NECESSARY 
REM helper programs LATER.EXE and SETERR.EXE are used

REM Check if EXE, AXS and AXT are later than TAB and STI; if so, skip
LATER %1.tab %1.sti / %1.exe %1.axt %1.axs
if errorlevel 1 goto skip
REM One of TAB or STI is later than EXE, AXS and AXT, so rerun TABLO and LTG
del %1.ax?
del %1.for
del %1.exe
echo on
tablo -sti %1.sti  >tb%1.log

if errorlevel 1 goto error
call ltg %1
if errorlevel 1 goto error
dir %1.exe
echo SUCCESSFULLY COMPILED %1
echo off
REM clean up junk files
del *.for
del *.lib
del *.mod
del modtable.txt
del opt
del opt90
del opt95
seterr 0
goto endbat
:error
echo off
echo ###### ERROR: FAILED TO COMPILE %1 #####
echo Check log file below
dir tb%1.log
dir %1.inf
echo Please press CTRL-C to terminate batch job
pause
goto endbat
:skip
seterr 0
echo off
echo COMPILE IS NOT NEEDED: %1.exe is later than %1.TAB
dir/od %1.*
:endbat
echo on
