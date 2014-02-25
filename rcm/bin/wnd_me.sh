matr=${1:-"l3_10"}

input_dir=../input2/for_wnd
output_dir=../output/for_wnd
alg=${alg:-wnd}

# gdb --args \
./rcm.out \
    ${input_dir}/${matr}.csr  \
    ${output_dir}/${matr}     \
    --threshold 0.1           \
    --info_file ${output_dir}/${matr}_${alg}.log  \
    --png ${output_dir}/${matr}            \
    --original_only \
    --algorithm ${alg}
