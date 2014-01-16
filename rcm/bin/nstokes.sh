#!/bin/bash
# RCM v. 0.1 runall-tool

OUTDIR=output3

cd `dirname $0`
mkdir -p ../${OUTDIR}/

# Stokes matricies
    EPS=0.1
    for simp_matrix in `cat ../input2/_stokes.list`; do 
        echo ${simp_matrix}
        ./rcm.out ../input2/${simp_matrix}.csr          \
                  ../${OUTDIR}/${simp_matrix}           \
            --info_file ../${OUTDIR}/${simp_matrix}.log \
            --threshold $EPS                            \
            --graph_threshold $EPS                      \
            --png ../${OUTDIR}/${simp_matrix};
    done
fi

# Stokes BIG matricies
if [ $STOKESBIG -eq 1 ]; then
    EPS=0.1
    for simp_matrix in `cat ../input2/_stokes_big.list`; do 
#        echo ${simp_matrix}
        ./rcm.out ../input2/${simp_matrix}.csr          \
                  ../${OUTDIR}/${simp_matrix}           \
            --info_file ../${OUTDIR}/${simp_matrix}.log \
            --threshold $EPS                            \
            --graph_threshold $EPS                      \
            --png ../${OUTDIR}/${simp_matrix} &
    done
fi

