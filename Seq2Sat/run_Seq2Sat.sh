for i in *_R1_001.fastq.gz; do j=$(basename ${i} _R1_001.fastq.gz); ../Seq2Sat/seq2sat -i ${i} -I ${j}_R2_001.fastq.gz --loc snp_loc_untrimmed.txt --var snp --sex sexLoc.txt --prefix ${j} --htJetter 15 --hmPerH 90 --hmPerL 85; done;

