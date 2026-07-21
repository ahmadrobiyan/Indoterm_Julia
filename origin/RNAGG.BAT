call deljunk.bat
echo on
del aggmod.har
agghar premod.har aggmod.har aggsupp.har -p >agghar.log
if errorlevel 1 goto error
dir aggmod.har
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





