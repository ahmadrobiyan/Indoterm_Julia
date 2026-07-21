call deljunk.bat

REM   delete existing files that this program will produce
del dumsupp.har
del dginit.har
del dumsec.agg
del dumreg.agg
del init.har

call compile init
if errorlevel 1 goto error
echo on

init -cmf init.cmf
if errorlevel 1 goto error
if not exist regsupp.har copy dumsupp.har regsupp.har
if not exist reg.agg copy dumreg.agg reg.agg
if not exist sec.agg copy dumsec.agg sec.agg

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
