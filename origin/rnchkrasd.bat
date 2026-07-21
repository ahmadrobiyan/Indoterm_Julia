call deljunk.bat
del dgchkrasD.HAR;

call compile chkras
if errorlevel 1 goto error
echo on

chkras -cmf chkrasD.cmf
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
