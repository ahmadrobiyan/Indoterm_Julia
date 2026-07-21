del *.log
del *.bak
del thismap.*
del term.sl4
del aggmod.har
del aggsets.har
rem Aggregate TERM database
rem  Takes 2 parameters: 
rem   (1) name of sector agg file without suffix
rem   (2) name of region agg file without suffix

if "%3"=="" goto OK1   
echo Too many parameters: DoAgg requires just 2
goto error
:OK1
if not "%2"=="" goto OK2   
echo Too few parameters: DoAgg requires 2
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

del aggmod.har
del aggsets.har

del sec.agg
del reg.agg
copy %1.agg sec.agg
copy %2.agg reg.agg

call rnPremod
call rnChkmdA
call rnPreagg
md mapstore
copy thismap.* mapstore\%2.*
call rnAgg
call rnAggset
call rnChkmdB
call rnTerm
grep -i "Total cpu" term.log
call rnChkmdC
dir/od *.sl4
echo   %1   %2

echo Finished OK
dir agg*.har/od
goto endbat
:error
echo ERROR
echo DoAgg requires 2 parameters:  doagg sec12 reg12
:endbat

