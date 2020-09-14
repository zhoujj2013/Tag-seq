gunzip data/*.fq.gz
perl ../bin/run_guideseq.pl config.TEST.txt all > config.TEST.log 2>config.TEST.err 
