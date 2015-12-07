/*=========================================================
 * eig_omp.cpp - Diagonalises a stack of matrices using
 *               Lapack calls and OpenMP.
 *
 * [V,E]=eig_omp_omp(H);
 *     where H is n x n x l, 
 *           V is n x n x l 
 *       and E is n x l
 *
 * This is a MEX-file for MATLAB.
 * 
 * Original Author: M. D. Le  [duc.le@stfc.ac.uk]
 * $Revision$ ($Date$)
 *=======================================================*/

#if !defined(_WIN32)
#define dsyevr dsyevr_
#define zheevr zheevr_
#endif

#include <cfloat>
#include <cstring>
#include <cmath>
#include "mex.h"
#include "matrix.h"
#include "lapack.h"

#ifndef _OPENMP
void omp_set_num_threads(int nThreads) {};
#define omp_get_num_threads() 1
#define omp_get_max_threads() 1
#define omp_get_thread_num()  0
#else
#include <omp.h>
#endif

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    mwSignedIndex m, n, nd;
    const size_t *dims;
    int ib, nblock, nb;
    int *blkid;
    char jobz;
    bool *issym, anynonsym=false;
    // Tolerance on whether matrix is symmetric/hermitian
    double tolsymm = sqrt(DBL_EPSILON);
    int nthread = omp_get_max_threads();
//  mexPrintf("Number of threads = %d\n",nthread);

    // Checks inputs
    if(nrhs!=1) {
        mexErrMsgIdAndTxt("eig_omp:nargin","Number of input argument must be 1.");
    }
    if(!mxIsNumeric(prhs[0])) {
        mexErrMsgIdAndTxt("eig_omp:notnumeric","Input must be a numeric array.");
    }
    nd = mxGetNumberOfDimensions(prhs[0]);
    if(nd<2 || nd>3) {
        mexErrMsgIdAndTxt("eig_omp:dims","Only 2D or 3D arrays are supported.");
    }
    dims = mxGetDimensions(prhs[0]);
    m = dims[0];
    n = dims[1];
    if(m!=n) {
        mexErrMsgIdAndTxt("eig_omp:square","Input matrix is not square.");
    }
    if(nd==3)
        nblock = (int)dims[2];
    else
        nblock = 1;

    // More efficient to group blocks together to run in a single thread than to spawn one thread per matrix.
    if(nblock<nthread*10) {
        nthread = 1;
    }
    blkid = new int[nthread+1];
    blkid[0] = 0;
    blkid[nthread] = nblock;
    nb = nblock/nthread;
    for(int nt=1; nt<nthread; nt++) {
       blkid[nt] = blkid[nt-1]+nb;
    }

    // Checks if all matrices are symmetric / hermitian.
    issym = new bool[nblock];
    // Initially assume all matrices symmetric 
    memset(issym,true,nblock*sizeof(bool));

#pragma omp parallel default(none) shared(issym,prhs) firstprivate(nthread, m, tolsymm, ib, blkid)
    if(mxIsComplex(prhs[0])) {
        double *A = mxGetPr(prhs[0]);
        double *Ai = mxGetPi(prhs[0]);
#pragma omp for
        for(int nt=0; nt<nthread; nt++) {
            for(ib=blkid[nt]; ib<blkid[nt+1]; ib++) {
                for(int ii=0; ii<m; ii++) {
                    // Just need to iterate over upper triangle
                    for(int jj=ii; jj<m; jj++) {
                        if(fabs( *(A+ii+jj*m) - *(A+jj+ii*m) ) > tolsymm
                         ||fabs( *(Ai+ii+jj*m) + *(Ai+jj+ii*m) ) > tolsymm) {
                            issym[ib] = false;
                            break;
                        }
                    }
                    if(!issym[ib])
                        break;
                }
                A += m*m;
                Ai += m*m;
            }
        }
    }
    else {
        double *A = mxGetPr(prhs[0]);
#pragma omp for
        for(int nt=0; nt<nthread; nt++) {
            for(ib=blkid[nt]; ib<blkid[nt+1]; ib++) {
                for(int ii=0; ii<m; ii++) {
                    for(int jj=ii; jj<m; jj++) {
                        if(fabs( *(A+ii+jj*m) - *(A+jj+ii*m) ) > tolsymm) {
                            issym[ib] = false;
                            break;
                        }
                    }
                    if(!issym[ib])
                        break;
                }
                A += m*m;
            }
        }
    }
    for(ib=0; ib<nblock; ib++)
        if(!issym[ib]) {
            anynonsym = true;
            break;
        }

    // Creates outputs
    if(nlhs<=1) {
        jobz = 'N';
        // Some matrices are not symmetric, use [ZD]GEEV algorithm, get complex eigenvalues
        if(anynonsym) {    
            if(nd==2)
                plhs[0] = mxCreateDoubleMatrix(m, 1, mxCOMPLEX);
            else
                plhs[0] = mxCreateDoubleMatrix(m, nblock, mxCOMPLEX);
        }
        else {
            if(nd==2)
                plhs[0] = mxCreateDoubleMatrix(m, 1, mxREAL);
            else
                plhs[0] = mxCreateDoubleMatrix(m, nblock, mxREAL);
        }
    }
    else {
        jobz = 'V';
        if(anynonsym) {    
            if(nd==2)
                plhs[1] = mxCreateDoubleMatrix(m, m, mxCOMPLEX);
            else
                plhs[1] = mxCreateDoubleMatrix(m, nblock, mxCOMPLEX);
        }
        else {
            if(nd==2)
                plhs[1] = mxCreateDoubleMatrix(m, m, mxREAL);
            else
                plhs[1] = mxCreateDoubleMatrix(m, nblock, mxREAL);
        }
    }
    // If some matrices are not symmetric, will get complex conjugate eigenpairs
    if(mxIsComplex(prhs[0]) || anynonsym) {
        if(nlhs>1) {
            if(nd==2)
                plhs[0] = mxCreateDoubleMatrix(m, m, mxCOMPLEX);
            else
                plhs[0] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxCOMPLEX);
        }
    }
    else {
        if(nlhs>1) {
            if(nd==2)
                plhs[0] = mxCreateDoubleMatrix(m, m, mxREAL);
            else
                plhs[0] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
        }
    }

#pragma omp parallel default(none) shared(plhs,prhs) firstprivate(nthread, m, nlhs, nd, ib, blkid, jobz, anynonsym, issym)
    {
#pragma omp for
        for(int nt=0; nt<nthread; nt++) {
            // These variables must be declared within the loop to make them local (and private)
            //   or we get memory errors.
            double *M, *E, *V, *D, *Di, *ptr_M, *ptr_Mi, *ptr_V, *ptr_Vi;
            mwSignedIndex m2, m22;
            size_t msz;
            char uplo = 'U';
            char range = 'A';
            char jobzn = 'N';
            double vl = -DBL_MAX, vu = DBL_MAX;
            mwSignedIndex il = 0, iu;
            double abstol = sqrt(DBL_EPSILON);
            mwSignedIndex lda, ldz, numfind;
            mwSignedIndex info, lwork, liwork, lzwork;
            mwSignedIndex *isuppz, *iwork;
            double *work, *zwork;
            int ii, jj;
            lda = m;
            ldz = m;
            iu = m;
            m2 = m*m;
            m22 = 2*m2;
            msz = m2*sizeof(double);
            lwork = 26*m;
            liwork = 10*m;
            isuppz = new mwSignedIndex[2*m];
            work = new double[lwork];
            iwork = new mwSignedIndex[liwork];
            if(mxIsComplex(prhs[0])) {
                lzwork = 4*m;
                M = new double[m22];
                zwork = new double[lzwork*2];
                if(nlhs>1)
                    V = new double[m22];
            }
            else
                M = new double[m2];
            // The output of the _evr Lapack routines gives eigenvalues as vectors
            // If we want it as a diagonal matrix, need to use a temporary array...
            if(mxIsComplex(prhs[0]) && anynonsym) {
                D = new double[2*m];
            }
            else if(nlhs>1 && nd==2) {
                D = new double[m];
                if(anynonsym)
                    Di = new double[m];
            }
            // Actual loop over individual matrices start here
            for(ib=blkid[nt]; ib<blkid[nt+1]; ib++) {
                if(!(mxIsComplex(prhs[0]) && anynonsym)) {
                    if(nlhs<=1)
                        D = mxGetPr(plhs[0]) + ib*m;
                    else if(nd==3)
                        D = mxGetPr(plhs[1]) + ib*m;
                }
                if(mxIsComplex(prhs[0])) {
                    ptr_M = mxGetPr(prhs[0]) + ib*m2;
                    ptr_Mi = mxGetPi(prhs[0]) + ib*m2;
                    // Interleaves complex matrices - Matlab stores complex matrix as an array of real
                    //   values followed by an array of imaginary values; Fortran (and C++ std::complex)
                    //   and hence Lapack stores it as arrays of pairs of values (real,imaginary).
                    for(ii=0; ii<m; ii++) {
                        for(jj=0; jj<m; jj++) {
                            M[ii*2+jj*2*m] = *ptr_M++;
                            M[ii*2+jj*2*m+1] = *ptr_Mi++;
                        }
                    }
                    // This does the actual computation - call using Matlab-included Intel MKL.
                    if(anynonsym)
                        zgeev(&jobz, &jobzn, &m, M, &lda, D, V, &ldz, V, &ldz, zwork, &lzwork, work, &info);
                    else
                        zheevr(&jobz, &range, &uplo, &m, M, &lda, &vl, &vu, &il, &iu, &abstol, &numfind, 
                                D, V, &ldz, isuppz, zwork, &lzwork, work, &lwork, iwork, &liwork, &info);
                    if(nlhs>1) {
                        ptr_V = mxGetPr(plhs[0]) + ib*m2;
                        ptr_Vi = mxGetPi(plhs[0]) + ib*m2;
                        for(ii=0; ii<m; ii++) {
                            for(jj=0; jj<m; jj++) {
                                *ptr_V++ = V[ii*2*m+jj*2];
                                *ptr_Vi++ = -V[ii*2*m+jj*2+1];
                            }
                        }
                    }
                }
                else {
                    ptr_M = mxGetPr(prhs[0]) + ib*m2;
                    V = mxGetPr(plhs[0]) + ib*m2;
                    if(anynonsym) {
                        if(nlhs<=1)
                            Di = mxGetPi(plhs[0]) + ib*m;
                        else if(nd==3)
                            Di = mxGetPi(plhs[1]) + ib*m;
                    }
                    // We cannot pass the right-side argument directly to Lapack because the routine
                    //   overwrites the input matrix to be diagonalised - so we create a copy
                    //   Matlab uses column-major like Fortran/Lapack (unlike C++) so can just memcpy
                    memcpy(M,ptr_M,msz);
                    // This does the actual computation - call using Matlab-included Intel MKL.
                    if(anynonsym)
                        dgeev(&jobzn, &jobz, &m, M, &lda, D, Di, V, &ldz, V, &ldz, work, &lwork, &info);
                    else
                        dsyevr(&jobz, &range, &uplo, &m, M, &lda, &vl, &vu, &il, &iu, &abstol, &numfind, 
                                D, V, &ldz, isuppz, work, &lwork, iwork, &liwork, &info);
                    // Need to account for complex conjugate eigenvalue/vector pairs (imaginary parts 
                    //   of eigenvectors stored in consecutive columns)
                    if(nlhs>1 && !issym[ib]) {
                        ptr_V = mxGetPr(plhs[0]) + ib*m2;
                        ptr_Vi = mxGetPi(plhs[0]) + ib*m2;
                        // Complex conjugate pairs of eigenvalues appear consecutively with the eigenvalue
                        //   having the positive imaginary part first.
                        for(ii=0; ii<m; ii++) {
                            if(*(Di+ii)>abstol) {
                                for(jj=0; jj<m; jj++) {
                                    *(ptr_Vi+(ii+1)*m+jj) = -*(ptr_V+(ii+1)*m+jj);
                                    *(ptr_Vi+ii*m+jj) = *(ptr_V+(ii+1)*m+jj);
                                    *(ptr_V+(ii+1)*m+jj) = *(ptr_V+ii*m+jj);
                                }
                            }
                        }
                    }
                }
                // If the output needs eigenvalues as diagonal matrix, creates it.
                if(nlhs>1 && nd==2) { 
                    if(mxIsComplex(prhs[0]) && anynonsym) {
                        E = mxGetPr(plhs[1]) + ib*m2;
                        for(ii=0; ii<m; ii++) 
                            *(E+ii+ii*m) = *(D+2*ii);
                        E = mxGetPi(plhs[1]) + ib*m2;
                        for(ii=0; ii<m; ii++) 
                            *(E+ii+ii*m) = *(D+2*ii+1);
                    }
                    else {
                        E = mxGetPr(plhs[1]) + ib*m2;
                        for(ii=0; ii<m; ii++) 
                            *(E+ii+ii*m) = *(D+ii);
                        if(anynonsym) {
                            E = mxGetPi(plhs[1]) + ib*m2; 
                            for(ii=0; ii<m; ii++) 
                                *(E+ii+ii*m) = *(Di+ii);
                        }
                    }
                }
                else if(mxIsComplex(prhs[0]) && anynonsym) {
                    if(nlhs<=1)
                        E = mxGetPr(plhs[0]) + ib*m;
                    else if(nd==3)
                        E = mxGetPr(plhs[1]) + ib*m;
                    for(ii=0; ii<m; ii++) 
                        *(E+ii) = *(D+2*ii);
                    if(nlhs<=1)
                        E = mxGetPi(plhs[0]) + ib*m;
                    else if(nd==3)
                        E = mxGetPi(plhs[1]) + ib*m;
                    for(ii=0; ii<m; ii++) 
                        *(E+ii) = *(D+2*ii+1);
                }
            }
            // Free memory...
            delete[]work; delete[]iwork; delete[]isuppz; delete[]M;
            if(mxIsComplex(prhs[0]))  {
                delete[]zwork;
                if(nlhs>1)
                    delete[]V;
            }
            if(nlhs>1 && nd==2)
                delete[]D;
        }
    }
    delete[]blkid; delete[]issym;
}
