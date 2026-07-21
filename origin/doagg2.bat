del *.log
del *.bak
del term.sl4
del thismap.*
del aggmod.har
del aggsets.har

rem Aggregate TERM database
rem  Takes 3 parameters: 
rem   (1) name of sector agg file without suffix
rem   (2) name of region agg file without suffix
rem   (3) name used for agg_%3.zip

if "%4"=="" goto OK1   
echo Too many parameters: DoAgg requires just 3
goto error
:OK1
if not "%3"=="" goto OK2   
echo Too few parameters: DoAgg requires 3
goto error
:OK2
if exist %1.agg goto OK3    
echo File %1.agg does not exist
goto error
:OK3
if exist %2.agg goto OK4    
echo File %2.agg does not exist
goto error
:OK4
del %3.zip 
call doagg %1 %2
zip agg_%3 aggmod.har aggsets.har
copy term.sl4 agg_%3.sl4
dir %3.zip 
echo Finished OK
goto endbat
:error
echo ERROR
echo DoAgg2 requires 3 parameters:  doagg sec12 reg12 zipname
:endbat

