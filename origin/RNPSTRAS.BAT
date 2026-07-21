call deljunk.bat
del chkrasall.har
del repeat.dat
echo on
call compile pstras
if errorlevel 1 goto error
echo on
copy distgone.har distgone.old
pstras -cmf pstras.cmf 
if errorlevel 1 goto error
cmbhar -sti cmbhar.sif
tsthar dgpstras.har flag
if errorlevel 2 goto error
if errorlevel 1 echo PLEASE REPEAT >repeat.dat
dir/od *.har
echo BATCH JOB SUCCESSFUL
dir pstras.har
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
