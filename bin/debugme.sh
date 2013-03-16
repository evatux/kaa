if [ -z "$indir" ]; then
    indir=nstokes
fi
if [ -z "$matr" ]; then
matr=nstokes_3_1
fi
if [ -z "$alg" ]; then
alg=rcm
fi

./rcm.out \
    ../input2/${indir}/${matr}.csr \
    ../output/${alg}/${matr}     \
    --threshold 0.1          \
    --graph_threshold 1      \
    --cheps_threshold 1e-6	\
    --info_file ../output/${alg}/${matr}.log  \
    --png ../output/${alg}/${matr}            \
    --algorithm ${alg}

