call deljunk.bat
call compile reg2
if errorlevel 1 goto error
echo on
reg2 -cmf reg2.cmf 
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
