bsub -b -I -q q_sw_expr -N 1 -np 1 -cgsp 64 -J lxy_test -share_size 14000 ./RabbitQCPlus_mpi -i ../SRR2496709_1.fastq -I ../SRR2496709_2.fastq -w 1 -o p1.fq -O p2.fq
