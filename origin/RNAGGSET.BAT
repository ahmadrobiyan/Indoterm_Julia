call deljunk.bat
call compile aggset
if errorlevel 1 goto error
echo on
aggset -cmf aggset.cmf 
if errorlevel 1 goto error
dir/od *.har
echo BATCH JOB SUCCESSFUL
dir/od *.har
goto endbat
:error
echo off
echo ###### ERROR: BATCH JOB FAILED #####
echo Check log file; most recent is listed last
dir/od *.log
echo Please press CTRL-C to terminate batch job
pause
:endbat
