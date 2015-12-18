/*===========================================================
 * chol_omp.cpp - Forms the Cholesky factorisation of a stack
 *                of matrices using Lapack calls and OpenMP.
 *
 * R = chol_omp(M);             - Errors if M != +ve def
 * [R,p] = chol_omp(M);         - No error, R is p x p
 * L = chol_omp(M,'lower');
 * [L,p] = chol_omp(M,'lower');
 *     where M is n x n x l, M=R'*R=L*L'
 *           
 * Unlike the matlab built-in this mex does not handle sparse
 * matrices. In addition, spinW specific option are:
 *
 * ... = chol_omp(...,'tol',val)
 *     where val is a constant to add to get +ve def-ness.
 * [K2,invK] = chol_omp(M,'Colpa')
 *     where K2 = (R*gComm*R') is Hermitian.
 *           gComm = [1...1,-1...1] is the commutator
 *           invK = inv(R)
 *
 * Note this function does not check for symmetry.
 * Also note that we do not truncate the output matrix if
 * the input is not positive definite, unlike the built-in.
 * i.e. R or L is _always_ of size (n x n) not (p x p)
 *
 * This is a MEX-file for MATLAB.
 *
 * Original Author: M. D. Le  [duc.le@stfc.ac.uk]
 * $Revision$ ($Date$)
 *=========================================================*/

#include <cfloat>
#include <cstring>
#include <cmath>
#include "mex.h"
#include "matrix.h"
#include "blas.h"
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
    int ib, nblock, nb, err_code=0;
    bool do_Colpa = false;
    int *blkid;
    char uplo = 'U', *parstr;
    int nthread = omp_get_max_threads();
    double tol = 0.;
//  mexPrintf("Number of threads = %d\n",nthread);

    // Checks inputs
    if(!mxIsNumeric(prhs[0])) {
        mexErrMsgIdAndTxt("chol_omp:notnumeric","Input matrix must be a numeric array.");
    }
    nd = mxGetNumberOfDimensions(prhs[0]);
    if(nd<2 || nd>3) {
        mexErrMsgIdAndTxt("chol_omp:wrongdims","Only 2D or 3D arrays are supported.");
    }
    dims = mxGetDimensions(prhs[0]);
    m = dims[0];
    n = dims[1];
    if(m!=n) {
        mexErrMsgIdAndTxt("chol_omp:notsquare","Input matrix is not square.");
    }
    if(nd==3)
        nblock = (int)dims[2];
    else
        nblock = 1;
    // Checks for optional arguments
    for(ib=1; ib<nrhs; ib++)
    {
        if(mxIsChar(prhs[ib])) {
            parstr = mxArrayToString(prhs[ib]);
            if(strcmp(parstr,"lower")==0)
                uplo = 'L';
            else if(strcmp(parstr,"Colpa")==0)
                do_Colpa = true;
            else if(strcmp(parstr,"tol")==0) {
                if((nrhs-1)>ib) {
                    tol = *(mxGetPr(prhs[ib+1]));
                    ib++;
                }
                else
                    mexErrMsgIdAndTxt("chol_omp:badtolarg","'tol' option requires a scalar tolerance value.");
            }
        }
    }

    // More efficient to group blocks together to run in a single thread than to spawn one thread per matrix.
    if(nblock<nthread) {
        nthread = nblock;
    }
    blkid = new int[nthread+1];
    blkid[0] = 0;
    blkid[nthread] = nblock;
    nb = nblock/nthread;
    for(int nt=1; nt<nthread; nt++) {
       blkid[nt] = blkid[nt-1]+nb;
    }

    // Creates outputs
    if(mxIsComplex(prhs[0])) {
        if(nd==2)
            plhs[0] = mxCreateDoubleMatrix(m, m, mxCOMPLEX);
        else
            plhs[0] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxCOMPLEX);
        if(nlhs>1) {
            if(do_Colpa) {
                if(nd==2)
                    plhs[1] = mxCreateDoubleMatrix(m, m, mxCOMPLEX);
                else
                    plhs[1] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxCOMPLEX);
            }
            else {
                if(nd==2)
                    plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
                else
                    plhs[1] = mxCreateDoubleMatrix(1, nblock, mxREAL);
            }
        }
    }
    else {
        if(nd==2)
            plhs[0] = mxCreateDoubleMatrix(m, m, mxREAL);
        else
            plhs[0] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
        if(nlhs>1) {
            if(do_Colpa) {
                if(nd==2)
                    plhs[1] = mxCreateDoubleMatrix(m, m, mxREAL);
                else
                    plhs[1] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
            }
            else {
                if(nd==2)
                    plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
                else
                    plhs[1] = mxCreateDoubleMatrix(1, nblock, mxREAL);
            }
        }
    }

#pragma omp parallel default(none) shared(plhs,prhs,err_code) \
    firstprivate(nthread, m, nlhs, nd, ib, blkid, uplo, tol, do_Colpa)
    {
#pragma omp for
        for(int nt=0; nt<nthread; nt++) {
            // These variables must be declared within the loop to make them local (and private)
            //   or we get memory errors.
            double *M, *Mp, *ptr_M, *ptr_Mi, *ptr_I, *gComm;
            mwSignedIndex lda = m, m2 = m*m, info, ii, jj, kk;
            mwSignedIndex f = mxIsComplex(prhs[0]) ? 2 : 1;
            char diag = 'N';
            char trans = 'C';
            char side = 'R';
            double *alpha;

            if(mxIsComplex(prhs[0]))
                M = new double[2*m*m];
            if(do_Colpa) {
                if(mxIsComplex(prhs[0])) {
                    Mp = new double[2*m*m];
                    alpha = new double[2];
                }
                else {
                    Mp = new double[m*m];
                    alpha = new double[1];
                }
                alpha[0] = 1.;
            }

            // Actual loop over individual matrices start here
            for(ib=blkid[nt]; ib<blkid[nt+1]; ib++) {
                // Loop is in case we want to try again with a constant added to the diagonal
                for(kk=0; kk<(tol>0?2:1); kk++) {
                    // Populate the matrix input array (which will be overwritten by the Lapack function)
                    if(mxIsComplex(prhs[0])) {
                        memset(M, 0, 2*m*m*sizeof(double));
                        ptr_M = mxGetPr(prhs[0]) + ib*m2;
                        ptr_Mi = mxGetPi(prhs[0]) + ib*m2;
                        // Interleaves complex matrices - Matlab stores complex matrix as an array of real
                        //   values followed by an array of imaginary values; Fortran (and C++ std::complex)
                        //   and hence Lapack stores it as arrays of pairs of values (real,imaginary).
                        for(ii=0; ii<m; ii++) {
                            for(jj=(uplo=='U'?0:ii); jj<(uplo=='U'?ii+1:m); jj++) {
                                M[ii*2*m+jj*2] = ptr_M[ii*m+jj];
                                M[ii*2*m+jj*2+1] = ptr_Mi[ii*m+jj];
                            }
                        }
                    }
                    else {
                        // *potrf overwrites the input array - copy only upper or lower triangle of input.
                        M = mxGetPr(plhs[0]) + ib*m2;
                        ptr_M = mxGetPr(prhs[0]) + ib*m2;
                        if(uplo=='U')
                            for(ii=0; ii<m; ii++) 
                                memcpy(M+ii*m, ptr_M+ii*m, (ii+1)*sizeof(double));
                        else
                            for(ii=0; ii<m; ii++) 
                                memcpy(M+ii*m+ii, ptr_M+ii*m+ii, (m-ii)*sizeof(double));
                    }
                    // Matrix not positive definite - try adding a small number to the diagonal.
                    if(kk==1)
                        for(ii=0; ii<m; ii++)
                            M[ii*m*f+ii*f] += tol; 
                    // Calls the actual Lapack algorithm for real or complex input.
                    if(mxIsComplex(prhs[0]))
                        zpotrf(&uplo, &m, M, &lda, &info);
                    else
                        dpotrf(&uplo, &m, M, &lda, &info);
                    // Matrix is positive definite, break out of loop.
                    if(info==0)
                        break;
                }
                if(info>0) {
                    if(nlhs<=1 || do_Colpa) {
                        #pragma omp critical 
                        {
                            err_code = 1;
                        }
                        break;
                    }
                    else {
                        ptr_I = mxGetPr(plhs[1]) + ib;
                        *ptr_I = (double)info;
                        // Zeros the non positive parts of the factor.
                      //kk = (mwSignedIndex)info-1;
                      //if(uplo=='U')
                      //    for(ii=kk; ii<m; ii++)
                      //        memset(M+ii*f*m, 0, mxIsComplex(prhs[0]) ? m*2*sizeof(double) : m*sizeof(double));
                      //else
                      //    for(ii=0; ii<m; ii++)
                      //        memset(M+ii*f*m+kk*f, 0, mxIsComplex(prhs[0]) ? (m-kk)*2*sizeof(double) : (m-kk)*sizeof(double));
                    }
                }
                if(do_Colpa) {
                    // Computes the Hermitian K^2 matrix = R*gComm*R';
                    memcpy(Mp, M, mxIsComplex(prhs[0]) ? m*m*2*sizeof(double) : m*m*sizeof(double));
                    // Applies the commutator [1..1,-1..-1] to the cholesky factor transposed
                    for(ii=m/2; ii<m; ii++)
                        for(jj=0; jj<m; jj++)
                            Mp[ii*f*m+jj*f] = -Mp[ii*f*m+jj*f];
                    if(mxIsComplex(prhs[0]))
                        for(ii=m/2; ii<m; ii++)
                            for(jj=0; jj<m; jj++)
                                Mp[ii*f*m+jj*f+1] = -Mp[ii*f*m+jj*f+1];
                    if(mxIsComplex(prhs[0]))
                        ztrmm(&side, &uplo, &trans, &diag, &m, &m, alpha, M, &lda, Mp, &lda);
                    else
                        dtrmm(&side, &uplo, &trans, &diag, &m, &m, alpha, M, &lda, Mp, &lda);
                    // Computes the inverse of the triangular factors (still stored in M)
                    if(nlhs>1) {
                        if(mxIsComplex(prhs[0])) {
                            ztrtri(&uplo, &diag, &m, M, &lda, &info);
                        }
                        else {
                            M = mxGetPr(plhs[1]) + ib*m2;
                            ptr_M = mxGetPr(plhs[0]) + ib*m2;
                            if(uplo=='U')
                                for(ii=0; ii<m; ii++) 
                                    memcpy(M+ii*m, ptr_M+ii*m, (ii+1)*sizeof(double));
                            else
                                for(ii=0; ii<m; ii++) 
                                    memcpy(M+ii*m+ii, ptr_M+ii*m+ii, (m-ii)*sizeof(double));
                            dtrtri(&uplo, &diag, &m, M, &lda, &info);
                        }
                        if(info>0) {
                            #pragma omp critical 
                            {
                                err_code = 2;
                            }
                            break;
                        }
                        if(mxIsComplex(prhs[0])) {
                            ptr_M = mxGetPr(plhs[1]) + ib*m2;
                            ptr_Mi = mxGetPi(plhs[1]) + ib*m2;
                            for(ii=0; ii<m; ii++) {
                                for(jj=(uplo=='U'?0:ii); jj<(uplo=='U'?ii+1:m); jj++) {
                                    ptr_M[ii*m+jj] = M[ii*2*m+jj*2];
                                    ptr_Mi[ii*m+jj] = M[ii*2*m+jj*2+1];
                                }
                            }
                        }
                    }
                    // Finally copies the output of the K^2=K*g*K' calculation to Matlab left-hand-side
                    //   (We had to keep it memory to compute the inverse)
                    if(mxIsComplex(prhs[0])) {
                        ptr_M = mxGetPr(plhs[0]) + ib*m2;
                        ptr_Mi = mxGetPi(plhs[0]) + ib*m2;
                        for(ii=0; ii<m; ii++) {
                            for(jj=0; jj<m; jj++) {
                                *ptr_M++ = Mp[ii*2*m+jj*2];
                                *ptr_Mi++ = Mp[ii*2*m+jj*2+1];
                            }
                        }
                    }
                    else {
                        memcpy(mxGetPr(plhs[0])+ib*m2, Mp, m*m*sizeof(double));
                    }
                }
                else if(mxIsComplex(prhs[0])) {
                    // Copy lapack output to matlab with interleaving
                    ptr_M = mxGetPr(plhs[0]) + ib*m2;
                    ptr_Mi = mxGetPi(plhs[0]) + ib*m2;
                    for(ii=0; ii<m; ii++) {
                      //for(jj=(uplo=='U'?ii:0); jj<(uplo=='U'?m:ii+1); jj++) { }
                        for(jj=0; jj<m; jj++) {
                            *ptr_M++ = M[ii*2*m+jj*2];
                            *ptr_Mi++ = M[ii*2*m+jj*2+1];
                        }
                    }
                }
                if(err_code!=0)
                    break;     // One of the threads got a singular or not pos def error - break loop here.
            }
            // Free memory...
            if(mxIsComplex(prhs[0]))
                delete[]M;
            if(do_Colpa) {
                delete[]Mp; delete[]alpha;
            }
            #ifndef _OPENMP
                if(err_code!=0)
                    break;
            #endif
        }
    }
    delete[]blkid;
    if(err_code==1) 
        mexErrMsgIdAndTxt("chol_omp:notposdef","The input matrix is not positive definite.");
    else if(err_code==2) 
        mexErrMsgIdAndTxt("chol_omp:singular","The input matrix is singular.");
}
