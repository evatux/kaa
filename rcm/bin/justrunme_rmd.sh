
output_dir=output_rmd
sz=4

mkdir -p ../${output_dir}/rcm
mkdir -p ../${output_dir}/md
mkdir -p ../${output_dir}/nd
mkdir -p ../${output_dir}/rmd

alg=rcm
./rcm.out \
    ../input2/nstokes/nstokes_${sz}_1.csr    \
    ../${output_dir}/${alg}/nstokes_${sz}_1  \
    --threshold 0                            \
    --cheps_threshold 1e-9                   \
    --cheps_substitute 1e-9                  \
    --info_file ../${output_dir}/${alg}/nstokes_${sz}_1.log  \
    --png ../${output_dir}/${alg}/nstokes_${sz}_1            \
    --algorithm ${alg} -o

alg=md
./rcm.out \
    ../input2/nstokes/nstokes_${sz}_1.csr    \
    ../${output_dir}/${alg}/nstokes_${sz}_1  \
    --threshold 0                            \
    --cheps_threshold 1e-9                   \
    --cheps_substitute 1e-9                  \
    --info_file ../${output_dir}/${alg}/nstokes_${sz}_1.log  \
    --png ../${output_dir}/${alg}/nstokes_${sz}_1            \
    --algorithm ${alg} -o

alg=nd
./rcm.out \
    ../input2/nstokes/nstokes_${sz}_1.csr    \
    ../${output_dir}/${alg}/nstokes_${sz}_1  \
    --threshold 0                            \
    --cheps_threshold 1e-9                   \
    --cheps_substitute 1e-9                  \
    --info_file ../${output_dir}/${alg}/nstokes_${sz}_1.log  \
    --png ../${output_dir}/${alg}/nstokes_${sz}_1            \
    --algorithm ${alg} -o


alg=rmd

for level in `seq 1 12`; do
    ./rcm.out \
        ../input2/nstokes/nstokes_${sz}_1.csr              \
        ../${output_dir}/${alg}/nstokes_${sz}_1_l${level}  \
        --threshold ${level}                               \
        --cheps_threshold 1e-9                             \
        --cheps_substitute 1e-9                            \
        --info_file ../${output_dir}/${alg}/nstokes_${sz}_1_l${level}.log  \
        --png ../${output_dir}/${alg}/nstokes_${sz}_1_l${level}            \
        --algorithm ${alg} -m
done
