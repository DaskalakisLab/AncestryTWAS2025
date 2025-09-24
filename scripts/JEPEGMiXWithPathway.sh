# -S /bin/bash
#$ -e /data/ajajoo/PEC/PECJPEGMIX/ResultsAllcisWithPathway202306//SCZ3_PGC.e
#$ -o /data/ajajoo/PEC/PECJPEGMIX/ResultsAllcisWithPathway202306//SCZ3_PGC.o
# -N SCZ3_PGC
# cd to the directory from which I submitted the
# job.  Otherwise it will execute in my data directory.
echo Beginning to run R job
date
/data/cchatzinakos/anaconda3/envs/r_env/bin/R CMD BATCH --no-restore --no-save --no-readline -iseed=1 -fr.aa=/data/ajajoo/PEC/PECJPEGMIX/ResultsAllcisWithPathway202306//SCZ3_PGC.genes.txt -fr.gg=/data/ajajoo/PEC/PECJPEGMIX/ResultsAllcisWithPathway202306//SCZ3_PGC.corr.txt -fr.bb=/data/ajajoo/PEC/PECJPEGMIX/ResultsAllcisWithPathway202306//SCZ3_PGC.agn.txt -fr.cc=/data/ajajoo/ASHG_TWAS/GWAS/SCZ/JEPEG/SCZ3_PGC.txt -fr.dd=/data/ajajoo/PEC/PECJPEGMIX/ResultsAllcisWithPathway202306//SCZ3_PGC.dir.txt -fr.ee=/data/ajajoo/PEC/PECJPEGMIX/ResultsAllcisWithPathway202306//logs/tpath_SCZ3_PGC.txt -fr.ss=/data/ajajoo/PEC/PECJPEGMIX/ResultsAllcisWithPathway202306//logs/ttpath_SCZ3_PGC.txt -nsim=2 -out=/data/ajajoo/PEC/PECJPEGMIX/ResultsAllcisWithPathway202306//logs/results_SCZ3_PGC.txt /data/ajajoo/PEC/PECJPEGMIX/script/20230628/callingJpegAllcisWithPathway202306.txt &> /data/ajajoo/PEC/PECJPEGMIX/ResultsAllcisWithPathway202306//SCZ3_PGC.Rout
date
sleep 120

