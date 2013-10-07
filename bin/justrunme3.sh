
alg=rcm
matr=bcsstm34
threshold=2e-5
cheps_threshold=1e-4
cheps_substitute=1e-4

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
