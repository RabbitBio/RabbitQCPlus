bsub -b -I -q q_sw_expr -n 1 -cgsp 64 -J lxy_test -share_size 7000 ./RabbitQCPlus -i ../r1.fq  -w 1 -o p1.fq
