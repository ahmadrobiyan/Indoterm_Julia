call deljunk.bat
del reg1.har
del diag1.har 


REM next lines implement bootstrap for average distance travelled
if not exist distgone.har copy distinit.har distgone.har

call compile reg1
if errorlevel 1 goto error
echo on

reg1 -cmf reg1.cmf
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
