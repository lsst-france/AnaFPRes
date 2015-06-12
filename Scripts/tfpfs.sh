#!/usr/local/bin/bash
C='X'
if [[ $# -gt 0 ]] ; then
  C=$1 ;
fi
echo '-------- tfpfs.sh C=' $C ' JOB_ID=' $JOB_ID ' JOB_NAME=' $JOB_NAME ' ---- '
## Initializing SOPHYA (to run the test programs)
source /sps/lsst/Library/Sophya/env.sh sl6 
BD=/sps/lsst/data/dev/lsstprod/DC_2013/forcedPhot_dir/forcedPhot/
EXED=/sps/lsst/dev/ansari/Ana_DC_LSST/Objs/
INP=/sps/lsst/dev/ansari/TstFPS/
OUTP=/sps/lsst/dev/ansari/TstFPS/
FILT='g'
${EXED}/fpfscan check -basedir $BD -outp ${OUTP}fpsumm_$C -camcols 1,2,3,4,5,6 -filters $FILT -runlistfile ${IN
P}runf_list_${C}.txt -fields 1,999 -prtlev 1 | tee $OUTP/tfpfs_${C}.log 
# ${EXED}/fpfscan check -basedir $BD -outp $OUTP -camcols 1 -filters i -runlistfile ${INP}micro_runf_list.txt -
fields 1,999 -prtlev 1 | tee $OUTP/tfpfs.log 

