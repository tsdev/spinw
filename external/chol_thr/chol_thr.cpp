/*===========================================================
 * chol_thr.cpp - Forms the Cholesky factorisation of a stack
 *                of matrices using Lapack calls and threads.
 *
 * R = chol_thr(M);             - Errors if M != +ve def
 * [R,p] = chol_thr(M);         - No error, R is p x p
 * L = chol_thr(M,'lower');
 * [L,p] = chol_thr(M,'lower');
 *     where M is n x n x l, M=R'*R=L*L'
 *           
 * Unlike the matlab built-in this mex does not handle sparse
 * matrices. In addition, spinW specific option are:
 *
 * ... = chol_thr(...,'tol',val)
 *     where val is a constant to add to get +ve def-ness.
 * [K2,invK] = chol_thr(M,'Colpa')
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

#define _THREADS
#if defined(__linux__)
#include <sys/sysinfo.h>
#elif defined(__FreeBSD__) || defined(__APPLE__)
#include <sys/types.h>
#include <sys/sysctl.h>
#else
#include <windows.h>
#endif

// We have to use threads directly rather than OpenMP because the Matlab included
//   Blas/Lapack libraries are compiles with OpenMP and are not compatible
#ifdef _WIN32
    #include <windows.h>
    #define MUTEX_LOCK     EnterCriticalSection
    #define MUTEX_UNLOCK   LeaveCriticalSection
    #define MUTEX_TYPE     CRITICAL_SECTION
    #define MUTEX_INIT(m)  InitializeCriticalSection (&m)
    #define EVENT_TYPE     HANDLE
    #define EVENT_INIT(e)  e = CreateEvent (NULL, TRUE, FALSE, NULL)
    #define EVENT_SIG(e)   SetEvent(e)
    #define THRLC_TYPE     DWORD
    #define THRLC_INIT(k)  k = TlsAlloc()
    #define THRLC_FREE(k)  TlsFree(k)
    #define THRLC_SET(k,v) TlsSetValue (k,v)
    #define THRLC_GET(v)   TlsGetValue (v)
    #define THRLC_GET_FAIL 0
#else
    #include <pthread.h>
    #define MUTEX_LOCK     pthread_mutex_lock
    #define MUTEX_UNLOCK   pthread_mutex_unlock
    #define MUTEX_TYPE     pthread_mutex_t
    #define MUTEX_INIT(m)  pthread_mutex_init (&m, NULL)
    #define EVENT_TYPE     pthread_cond_t
    #define EVENT_INIT(e)  pthread_cond_init (&e, NULL)
    #define EVENT_SIG(e)   pthread_cond_signal (&e)
    #define THRLC_TYPE     pthread_key_t
    #define THRLC_INIT(k)  pthread_key_create(&k, dataDestructor)
    #define THRLC_FREE(k)  pthread_key_delete(k)
    #define THRLC_SET(k,v) pthread_setspecific (k,v)
    #define THRLC_GET(v)   pthread_getspecific (v)
    #define THRLC_GET_FAIL NULL
    void dataDestructor(void *data) { }
#endif

// ----------------------------------------------------------------------------------- //
// Define a struct containing all required inputs to the threads which do not change.
// ----------------------------------------------------------------------------------- //
typedef struct {
    mwSignedIndex m;
    int nlhs;
    mwSignedIndex nd;
    int *blkid;
    char uplo;
    double tol;
    bool do_Colpa;
    mxArray **plhs;
    const mxArray **prhs;
    int *err_code;
    int *warn1;
} global_thread_data;
// ----------------------------------------------------------------------------------- //
// Class for inputs which differ from thread to thread.
// ----------------------------------------------------------------------------------- //
class thread_input {
    public:
        int nt;
        thread_input(int _nt)
        {
            nt = _nt;
        }
};

// ----------------------------------------------------------------------------------- //
// Declares these variables global, so all threads can see them
// ----------------------------------------------------------------------------------- //
global_thread_data gtd;
#define NUM_THREADS 64          // Can support up to a maximum of 64 threads
thread_input *tin[NUM_THREADS]; //   - hard coded because global variable.
MUTEX_TYPE mutex_error;
MUTEX_TYPE mutex_warn;
EVENT_TYPE checkfinish;
THRLC_TYPE threadSpecificKey;

//#pragma omp parallel default(none) shared(plhs,prhs,err_code,warn1) \
//    firstprivate(nthread, m, nlhs, nd, ib, blkid, uplo, tol, do_Colpa)
#ifdef _WIN32
DWORD WINAPI thread_iteration(void *input)
#else
void *thread_iteration(void *input)
#endif
{
    thread_input *td;
    mwSignedIndex m = gtd.m;
    int nlhs = gtd.nlhs;
    mwSignedIndex nd = gtd.nd;
    int *blkid = gtd.blkid;
    char uplo = gtd.uplo;
    double tol = gtd.tol;
    bool do_Colpa = gtd.do_Colpa;
    mxArray **plhs = gtd.plhs;
    const mxArray **prhs = gtd.prhs;
    int *err_code = gtd.err_code;
    int *warn1 = gtd.warn1;
    td = (thread_input*)input;
    int nt = td->nt;
    // These variables must be declared within the loop to make them local (and private)
    //   or we get memory errors.
    double *M, *Mp, *ptr_M, *ptr_Mi, *ptr_I, *gComm;
    mwSignedIndex lda = m, m2 = m*m, info, ii, jj, kk;
    mwSignedIndex f = mxIsComplex(prhs[0]) ? 2 : 1;
    char diag = 'N';
    char trans = 'C';
    char side = 'R';
    double *alpha;
    int ib;

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
            if(kk==1) {
              //#pragma omp critical 
                MUTEX_LOCK(&mutex_warn);
                {
                    *warn1 = 1;
                }
                MUTEX_UNLOCK(&mutex_warn);
                for(ii=0; ii<m; ii++)
                    M[ii*m*f+ii*f] += tol; 
            }
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
              //#pragma omp critical 
                MUTEX_LOCK(&mutex_error);
                {
                    *err_code = 1;
                }
                MUTEX_UNLOCK(&mutex_error);
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
                  //#pragma omp critical 
                    MUTEX_LOCK(&mutex_error);
                    {
                        *err_code = 2;
                    }
                    MUTEX_UNLOCK(&mutex_error);
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
        if(*err_code!=0)
            break;     // One of the threads got a singular or not pos def error - break loop here.
    }
    // Free memory...
    if(mxIsComplex(prhs[0]))
        delete[]M;
    if(do_Colpa) {
        delete[]Mp; delete[]alpha;
    }
    #if !defined(_THREADS) || defined(_WIN32)
    return 0;
    #endif
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    mwSignedIndex m, n, nd;
    const size_t *dims;
    int ib, nblock, nb, nt, err_code=0;
    bool do_Colpa = false;
    int *blkid;
    char uplo = 'U', *parstr;
    double tol = 0.;
    int warn1 = 0;
    int nthread;

    // Checks if user specified number of threads in a global variable
    mxArray *thrptr = mexGetVariable("global","sw_num_threads");

    // System-dependent calls to find number of processors (from GotoBLAS)
    if(thrptr == NULL)
    {
    #if defined(__linux__)
        nthread = get_nprocs();
    #elif defined(__FreeBSD__) || defined(__APPLE__)
        int m[2], count; size_t len;
        m[0] = CTL_HW; m[1] = HW_NCPU; len = sizeof(int);
        sysctl(m, 2, &nthread, &len, NULL, 0);
    #else
        SYSTEM_INFO sysinfo; GetSystemInfo(&sysinfo);
        nthread = sysinfo.dwNumberOfProcessors;
    #endif
    }
    else {
        nthread = (int)mxGetPr(thrptr)[0];
        mxDestroyArray(thrptr);
    }
    if(nthread>NUM_THREADS) nthread=NUM_THREADS;
    mexPrintf("Using %d threads.\n",nthread);

    // Checks inputs
    if(!mxIsNumeric(prhs[0])) {
        mexErrMsgIdAndTxt("chol_thr:notnumeric","Input matrix must be a numeric array.");
    }
    nd = mxGetNumberOfDimensions(prhs[0]);
    if(nd<2 || nd>3) {
        mexErrMsgIdAndTxt("chol_thr:wrongdims","Only 2D or 3D arrays are supported.");
    }
    dims = mxGetDimensions(prhs[0]);
    m = dims[0];
    n = dims[1];
    if(m!=n) {
        mexErrMsgIdAndTxt("chol_thr:notsquare","Input matrix is not square.");
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
                    mexErrMsgIdAndTxt("chol_thr:badtolarg","'tol' option requires a scalar tolerance value.");
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
    for(nt=1; nt<nthread; nt++) {
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

    // Sets values of the shared variables, accessible through the global struct gtd
    gtd.m = m;
    gtd.nlhs = nlhs;
    gtd.nd = nd;
    gtd.blkid = blkid;  //
    gtd.uplo = uplo;
    gtd.tol = tol;
    gtd.do_Colpa = do_Colpa;
    gtd.plhs = plhs;    //
    gtd.prhs = prhs;    //
    gtd.err_code = &err_code;  //
    gtd.warn1 = &warn1;  //
    MUTEX_INIT(mutex_error);
    MUTEX_INIT(mutex_warn);
    EVENT_INIT(checkfinish);
    THRLC_INIT(threadSpecificKey);
    #if defined  (__linux__) || defined (__APPLE__)
        pthread_t threads[NUM_THREADS]; int rc; void *status;
        pthread_attr_t attr;
        pthread_attr_init(&attr);
        pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
    #else
        HANDLE threads[NUM_THREADS];
        DWORD tid[NUM_THREADS], dwError;
    #endif
    // Starts threads or run through blocks of matrices in serial
    for(nt=0; nt<nthread; nt++) {
        tin[nt] = new thread_input(nt);
    #if defined  (__linux__) || defined (__APPLE__)
        rc = pthread_create(&threads[nt], &attr, thread_iteration, (void *)tin[nt]);
        if(rc) {
            mexErrMsgIdAndTxt("chol_thr:threadcreate","Cannot create thread %i due to error code %i",nt+1,rc);
        }
    #else
        threads[nt] = CreateThread(NULL, 0, thread_iteration, (void *)tin[nt], 0, &tid[nt]);
        if(threads[nt]==NULL) {
            dwError=GetLastError();
            mexErrMsgIdAndTxt("chol_thr:threadcreate","Cannot create thread %i due to error code %i",nt+1,dwError);
        }
    #endif
    }
    // Wait for all threads to finish and clean up
    for(nt=0; nt<nthread; nt++) {
    #if defined  (__linux__) || defined (__APPLE__)
        rc = pthread_join(threads[nt], &status);
        if(rc) {
            mexErrMsgIdAndTxt("chol_thr:threadjoin","Cannot end/join thread %i due to error code %i.",nt+1,rc);
        }
    #else
        if(WaitForSingleObject(threads[nt],INFINITE)==0xFFFFFFFF) {
            mexErrMsgIdAndTxt("chol_thr:threadjoin","Cannot end/join thread %i",nt+1);
        }
    #endif
        delete[]tin[nt];
    }
    #if defined (__linux__) || defined (__APPLE__)
    pthread_attr_destroy(&attr);
    pthread_mutex_destroy(&mutex_error);
    #endif
    THRLC_FREE(threadSpecificKey);

    delete[]blkid;
    if(err_code==1) 
        mexErrMsgIdAndTxt("chol_thr:notposdef","The input matrix is not positive definite.");
    else if(err_code==2) 
        mexErrMsgIdAndTxt("chol_thr:singular","The input matrix is singular.");
    if(warn1==1)
        mexWarnMsgIdAndTxt("chol_thr:NonPosDefHamiltonian","To make the Hamiltonian positive definite, a small omega_tol value was added to its diagonal!");
}
