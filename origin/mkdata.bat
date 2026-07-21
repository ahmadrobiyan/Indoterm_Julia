del *.log
del *.bak
del *.sl4

rem NB: lines beginning "rem" are comments only

rem activate [ie, un-comment-out] next line to do aggregation only
rem goto stage3

rem TERM data process
rem inputs are: National.HAR (compulsory)
rem and optionally:
rem   REGSUPP.HAR (additional regional data)
rem   DISTGONE.HAR (average origin-destination distance of commodities)  
rem   SEC.AGG REG.AGG (used by preagg.cmf)
rem All optional files are created (if not present) by the procedure

:STAGE1

call rnInit
call rnREG0

:STAGE2
 
call rnReg1
call rnReg2
call rnchkrasa
call rnTrdRAS
call rnchkrasb
rem call rnRASLIN
rem call rnchkrasc
call rnTrdRAS2
call rnchkrasd
call rnPstras

REM rnPstRas creates a file repeat.dat
REM  if the revised distance file distgone.har
REM  is much different from the previous file distgone.old.
REM If repeat.dat exists, STAGE2 is restarted
dir repeat.dat
if not exist repeat.dat goto norepeat
echo  Hit return to repeat average distance iteration
pause Average distance not yet converged
goto STAGE2
:norepeat
call rnPremod

:STAGE3
call rnChkmdA
rem rnPreAgg runs preagg.tab; also makes thismap.mif map for shademap
call rnPreagg
call rnAgg
call rnAggset
call rnChkmdB
call rnTerm
call rnChkmdC
dir/od *.sl4




