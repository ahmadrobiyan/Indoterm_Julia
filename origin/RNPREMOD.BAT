call deljunk.bat
call compile premod
if errorlevel 1 goto error
echo on
premod -cmf premod.cmf 
if errorlevel 1 goto error
dir/od *.har
echo BATCH JOB SUCCESSFUL
dir premod.har
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
