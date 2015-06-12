set LOGP = /sps/lsst/dev/ansari/TstFPS/
qsub -P P_lsst -l ct=4200,os=sl6,vmem=2G,sps=1 -o ${LOGP} -e ${LOGP}  ./tfpfs.sh
# qsub -P P_lsst -l ct=500,os=sl6,vmem=2G,sps=1 -o ${LOGP}tge_J2.log -e ${LOGP}tge_J2.log ./tgest.sh
qstat

