
alg=nd
matr=c-22
threshold=0.001
cheps_threshold=2e-3
cheps_substitute=2e-3

sz=3

./rcm.out \
    ../input2/fmc/${matr}.mtx                       \
    ../output4/${alg}/${matr}                       \
    --threshold ${threshold}                        \
    --cheps_threshold ${cheps_threshold}            \
    --cheps_substitute ${cheps_substitute}          \
    --info_file ../output4/${alg}/${matr}.log       \
    --png ../output4/${alg}/${matr}                 \
    --algorithm ${alg} --florida_file
