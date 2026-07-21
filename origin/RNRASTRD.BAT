call deljunk.bat
del rastrd.HAR;
del dgrastrd.HAR;

call compile RASTRD
if errorlevel 1 goto error
echo on



RASTRD -cmf RASTRD.cmf
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
