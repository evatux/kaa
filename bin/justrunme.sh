
alg=rcm
sz=3

./rcm.out \
    ../input2/nstokes/nstokes_${sz}_1.csr \
    ../output3/${alg}/nstokes_${sz}_1     \
    --threshold 0.1          \
    --cheps_threshold 1e-5          \
    --cheps_substitute 4.2e-3          \
    --info_file ../output3/${alg}/nstokes_${sz}_1.log  \
    --png ../output3/${alg}/nstokes_${sz}_1            \
    --algorithm ${alg}

