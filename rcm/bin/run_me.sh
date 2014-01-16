#!/bin/bash
# RCM v. 0.1 runall-tool

OUTDIR=output3

cd `dirname $0`
mkdir -p ../${OUTDIR}/

MATRIX_SET=1
STOKES=1
SIMPLE=1

while [ -n "$*" ]
do
    if [ "$1" == "stokes" ]; then
        if [ ${MATRIX_SET} -eq 1 ]; then
            MATRIX_SET=0;
            STOKES=1;
            SIMPLE=0;
            STOKESBIG=0;
        else
            STOKES=1;
        fi
    fi
    if [ "$1" == "stokes_big" ]; then
        if [ ${MATRIX_SET} -eq 1 ]; then
            MATRIX_SET=0;
            STOKES=0;
            SIMPLE=0;
            STOKESBIG=1;
        else
            STOKESBIG=1;
        fi
    fi
    if [ "$1" == "simple" ]; then
        if [ ${MATRIX_SET} -eq 1 ]; then
            MATRIX_SET=0;
            STOKES=0;
            SIMPLE=1;
            STOKESBIG=0;
        else
            SIMPLE=1;
        fi
    fi

    shift
done

# Simple matricies
if [ $SIMPLE -eq 1 ]; then
    EPS=1
    for simp_matrix in `cat ../input2/_simp.list`; do 
        ./rcm.out ../input2/${simp_matrix}.csr          \
                  ../${OUTDIR}/${simp_matrix}           \
            --info_file ../${OUTDIR}/${simp_matrix}.log \
            --threshold $EPS                            \
            --png ../${OUTDIR}/${simp_matrix};
    done
fi

# Stokes matricies
if [ $STOKES -eq 1 ]; then
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

