This download includes the full suite of data, programs and scripts used to generate the 2016 Indonesia TERM database
from a national CGE database (ORANI-G format) that includes the limited regional data needed for a "top-down" regional
extension. To replicate the data process, a source code GEMPACK licence is needed. Submitted JMH October 2021.

To run the data programs, unzip the archive into an empty new folder.
Open a CMD prompt in that folder and type:
MkData.bat 

You should see a sequence of programs running, concluding with a test simulation producing file Term.sl4.

The first time you run MkData.bat it will 
(1) run TABLO and LTG to make an EXE file for each TAB file.
(2) run the EXE files

If you run MkData.bat again, stage (1) will be omitted, saving time.

The starting point data files are:
NATIONAL.HAR          ORANI-G format with national data
REGSUPP.HAR           Additional regional detail 
sec.agg, reg.agg      HAR files specifying aggregation files 
DISTGONE.HAR..........Average origin-destination distance travelled by commodities (optional)

The file DatNotes.txt contains some notes on the programs used by mkdata.bat.
The file regdata.doc contains an overview of the process.

After you run MkData.bat, try running  DoAllAggs.bat: it produces six zips with different TERM aggregations, using DOAGG2.BAT and DOAGG.BAT
