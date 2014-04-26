SpMV framework

structure:
    common/         <-- common headers and core functions
    spmv/           <-- kernels
    tests/          <-- tests
        input/
        output/

API:
    int  spmv_desc_create (spmv_desc_t *desc, TMatrix_CSR *matr, void *info);
    int  spmv_compute     (spmv_desc_t *desc, TVector_SMP *in,   TVector_SMP *out);
    void spmv_desc_destroy(spmv_desc_t *desc);
