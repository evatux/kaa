
./rcm.out \
    ../input2/nstokes/nstokes_3_1.csr \
    ../output3/rcm/nstokes_3_1     \
    --threshold 0.1          \
    --cheps_threshold 1e-5          \
    --cheps_substitute 1e-5          \
    --info_file ../output3/rcm/nstokes_3_1.log  \
    --png ../output3/rcm/nstokes_3_1            \
    --algorithm rcm

