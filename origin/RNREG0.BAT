call deljunk.bat
REM   delete existing files that this program will produce
del reg0.har
del dgreg0.har
del natelast.har
del distinit.har 

call compile reg0
if errorlevel 1 goto error
echo on

reg0 -cmf reg0.cmf >nul:
if errorlevel 1 goto error

dir/od *.har
echo BATCH JOB SUCCESSFUL
goto endbat
:error
echo off
echo ###### ERROR: BATCH JOB FAILED #####
echo Check log file; most recent is listed last
dir/od *.log
:again
echo Please press CTRL-C to terminate batch job
pause > nul
goto again
:endbat

