if [ -z "$indir" ]; then
    indir=nstokes
fi
if [ -z "$matr" ]; then
matr=nstokes_5_2
fi
if [ -z "$alg" ]; then
alg=nd
fi

./rcm.out \
    ../input2/${indir}/${matr}.csr \
    ../output/${alg}/${matr}     \
    --threshold 0.1          \
    --graph_threshold 1      \
    --cheps_threshold 1e-6	\
    --cheps_substitute 0.0087 \
    --info_file ../output/${alg}/${matr}.log  \
    --png ../output/${alg}/${matr}            \
    --algorithm ${alg}

