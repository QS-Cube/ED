#PJM -L rscunit=bwmpc
#PJM -L rscgrp=batch
#PJM -L vnode=1
#PJM -L vnode-core=40
#PJM -L elapse=00:05:00
#PJM -g Q20458
#PJM -j
#
export OMP_NUM_THREADS=40
for i in {42..62}; do
echo $i
./QS3.exe < input/input_$i.dat > output/output_$i.dat
done
#
# pjsub --interact -L rscunit=bwmpc -L vnode=1 -L vnode-core=40 -g Q20458
