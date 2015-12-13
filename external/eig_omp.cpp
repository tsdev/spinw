/*===========================================================
 * eig_omp.cpp - Diagonalises a stack of matrices using
 *               Lapack calls and OpenMP.
 *
 * [V,E]=eig_omp(H,'orth','sort');
 *     where H is n x n x l, 
 *           V is n x n x l 
 *       and E is n x l
 *
 * 'orth' and 'sort' are options. See eig_omp.m for full info
 *
 * This is a MEX-file for MATLAB.
 *
 * In addition the the mexFunction, the following 
 * subsidiary functions are defined:
 *
 * void fliplr()    - Flips matrix by columns
 * void quicksort() - Quicksort algorithm, returns indices
 * void sort()      - Sort and permutes eigenvalue/vectors
 * void orth()      - Orthogonalises degenerate eigenvecs
 * 
 * Original Author: M. D. Le  [duc.le@stfc.ac.uk]
 * $Revision$ ($Date$)
 *=========================================================*/

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

// Flips a (column-major) matrix by columns, like the matlab function.
void fliplr(double *M, mwSignedIndex m, mwSignedIndex n, double *vec, bool isreal)
{
    int ii;
    double val, *p0r, *p1r, *p0i, *p1i;
    // Just a row vector, reverse order of values.
    if(m==1) {
        if(isreal) {
            p0r = M; p1r = M+n-1;
            for(ii=0; ii<(n/2); ii++) {
                val = *p1r;
                *p1r-- = *p0r;
                *p0r++ = val;
            }
        }
        else {
            p0r = M; p1r = M+n*2-2;
            p0i = M+1; p1i = M+n*2-1;
            for(ii=0; ii<(n/2); ii++) {
                val = *p1r; *p1r-- = *p0r; *p0r++ = val;
                val = *p1i; *p1i-- = *p0i; *p0i++ = val;
            }
        }
    }
    // Actual matrix - assume column major
    else {
        if(isreal) {
            size_t msz = m*sizeof(double);
            for(ii=0; ii<(n/2); ii++) {
                memcpy(vec, M+(n-ii-1)*n, msz);
                memcpy(M+(n-ii-1)*n, M+ii*n, msz);
                memcpy(M+ii*n, vec, msz);
            }
        }
        else {
            mwSignedIndex n2 = n*2;
            size_t m2sz = 2*m*sizeof(double);
            for(ii=0; ii<(n/2); ii++) {
                memcpy(vec, M+(n-ii-1)*n2, m2sz);
                memcpy(M+(n-ii-1)*n2, M+ii*n2, m2sz);
                memcpy(M+ii*n2, vec, m2sz);
            }
        }
    }
}

// Quicksort modified from public domain implementation by Darel Rex Finley.
//      http://alienryderflex.com/quicksort/
void quicksort(int *id, double *val, mwSignedIndex elements)
{
    double piv;
    int i=0, L, R, C, swap;
    int beg[300], end[300];
    beg[0]=0; end[0]=(int)elements;
    while (i>=0) {
        L=beg[i]; R=end[i]-1;
        if (L<R)
        {
            piv=val[id[L]]; C=id[L];
            while (L<R)
            {
                while (val[id[R]]>=piv && L<R) R--; if (L<R) id[L++]=id[R];
                while (val[id[L]]<=piv && L<R) L++; if (L<R) id[R--]=id[L];
            }
            id[L]=C; beg[i+1]=L+1; end[i+1]=end[i]; end[i++]=L;
            if (end[i]-beg[i]>end[i-1]-beg[i-1])
            {
                swap=beg[i]; beg[i]=beg[i-1]; beg[i-1]=swap;
                swap=end[i]; end[i]=end[i-1]; end[i-1]=swap;
            }
        }
        else i--;
    }
}

// This is called for real general matrices - complex conjugate eigenvectors stored in consecutive columns
void sort(mwSignedIndex m, double *Dr, double *Di, double *V, double *work, int sort_type)
{
    int *id, ii, jj;
    double *val, *pD, *pDi, *pV, abstol = sqrt(DBL_EPSILON);
    size_t msz;

    // Assume(!) that workspace is max[3*m+601,m*(m+3)] large.
    id = (int*)work;
    for(ii=0; ii<(int)m; ii++) 
        id[ii] = ii;

    // Probably should modify quicksort to allow descending sort
    if(sort_type==-1) {
        val = work+m+1;
        // For now just feed quicksort the negative of the eigenvalues...
        for(ii=0; ii<(int)m; ii++) 
            val[ii] = -Dr[ii];
    }
    else 
        val = Dr;

    // Does the sorting using QuickSort (probably not the quickest algorithm, but fast and simple...)
    quicksort(id, val, (int)m);

    // Now permute the eigenvalues and eigenvectors 
    //   this is horendously bloated... maybe try Fich et al. 1995 
    //   http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.29.2256
    msz = m*sizeof(double);
    pD = work+m; memcpy(pD, Dr, msz);
    pDi = work+2*m; memcpy(pDi, Di, msz);
    if(V!=0)
        pV = work+3*m; memcpy(pV, V, m*m*sizeof(double));
    for(ii=0; ii<(int)m; ii++) {
        // Take care to preserve the order of the eigenvectors (real parts first)
        if((ii+1)<(int)m)
            if(fabs(Dr[id[ii]]-Dr[id[ii+1]])<abstol && Di[id[ii]]<abstol) {
                jj = id[ii+1];
                id[ii+1] = id[ii];
                id[ii] = jj;
            }
        pD[ii] = Dr[id[ii]];
        pDi[ii] = Di[id[ii]];
        if(V!=0)
            memcpy(&pV[ii*m], &V[id[ii]*m], msz);
    }
    memcpy(Dr, pD, msz);
    memcpy(Di, pDi, msz);
    if(V!=0)
        memcpy(V, pV, m*m*sizeof(double));
}

// This is called for complex general matrices - both eigenvalues and eigenvectors are actually complex
void sort(mwSignedIndex m, double *D, double *V, double *work, int sort_type)
{
    int *id, ii;
    double *val, *pD, *pV;
    size_t msz;

    // Assume(!) that workspace is max[3*m+601,2*m*(m+3)] large.
    id = (int*)work;
    for(ii=0; ii<(int)m; ii++) 
        id[ii] = ii;

    // Populate the sort variable as the real part of the eigenvalues
    val = work+m+1;
    if(sort_type==-1) {
        for(ii=0; ii<(int)m; ii++) 
            val[ii] = -D[ii*2];
    }
    else {
        for(ii=0; ii<(int)m; ii++) 
            val[ii] = D[ii*2];
    }

    // Does the sorting using QuickSort
    quicksort(id, val, (int)m);

    // Now permute the eigenvalues and eigenvectors - in future look at using Fich et al. 1995
    msz = 2*m*sizeof(double);
    pD = work+m; memcpy(pD, D, msz);
    if(V!=0)
        pV = work+3*m; memcpy(pV, V, 2*m*m*sizeof(double));
    for(ii=0; ii<(int)m; ii++) {
        pD[ii*2] = D[id[ii]*2];
        pD[ii*2+1] = D[id[ii]*2+1];
        if(V!=0)
            memcpy(&pV[ii*2*m], &V[id[ii]*2*m], msz);
    }
    memcpy(D, pD, msz);
    if(V!=0)
        memcpy(V, pV, 2*m*m*sizeof(double));
}

// This is called for real general matrices - complex conjugate eigenvectors stored in consecutive columns
void orth(mwSignedIndex m, double *Dr, double *Di, double *Vr, double *Vi, double *work, bool isreal)
{
    int *id, ii, jj, kk, nn;
    mwSignedIndex n, info;
    double *pVz, *pS;
    double abstol = sqrt(DBL_EPSILON);
    char jobu = 'O';
    char jobvt = 'N';
    mwSignedIndex lwork = 3*m;
    double *zwork = work + (2*m*(m+7)) - 2*3*m;
    double *rwork = zwork - 5*m;

    // Assume workspace is size [m*(m+7)]*sizeof(complexdouble)
    id = (int*)work;
    for(ii=0; ii<(int)m; ii++) 
        id[ii] = ii;

    // Use QuickSort to get indices of eigenvalues in ascending order by real part
    quicksort(id, Dr, m);

    // Check for warning
//  for(ii=0; ii<(int)m; ii++)
//      if(fabs(Di[ii])>abstol) {
//          mexWarnMsgIdAndTxt("eig_omp:notdefinite","matrix contains complex eigenvalues, othogonalisation will not be accurate");
//          break;
//      }

    pVz = work+2*m;
    for(ii=0; ii<(int)m; ii++) {
        if((ii+1)<(int)m)
            if(fabs(Dr[id[ii]]-Dr[id[ii+1]])<abstol) {
                if(Di[id[ii]]<-abstol) {
                    jj = id[ii+1];
                    id[ii+1] = id[ii];
                    id[ii] = jj;
                }
                kk = ii;
                while(true)
                    if(fabs(Dr[id[ii]]-Dr[id[++kk]])>abstol) 
                        break;
                // Populates a complex matrix with the degenerate eigenvectors, taking care about interleaving
                n = 0;
                for(jj=ii; jj<kk; jj++) {
                    // Conjugate pairs appear consecutively with positive imaginary part first
                    if(Di[id[jj]]>abstol) {
                        for(nn=0; nn<m; nn++) {
                            pVz[n*2*m+nn*2] = Vr[id[jj]*m+nn];
                            pVz[n*2*m+nn*2+1] = Vr[id[jj+1]*m+nn];
                            pVz[(n+1)*2*m+nn*2] = Vr[id[jj]*m+nn];
                            pVz[(n+1)*2*m+nn*2+1] = -Vr[id[jj+1]*m+nn];
                        }
                        jj++;
                        n++;
                    }
                    else
                        for(nn=0; nn<m; nn++) 
                            pVz[n*2*m+nn*2] = Vr[id[jj]*m+nn];
                    n++;
                }
                // Does the singular value decomposition to get the orthogonal basis and singular values
                pS = work + 2*(n+1)*m;
                zgesvd(&jobu, &jobvt, &m, &n, pVz, &m, pS, pVz, &m, pVz, &m, zwork, &lwork, rwork, &info);
                // Checks that number of singular values == n
                for(nn=0; nn<n; nn++)
                    if(fabs(pS[nn])<0)
                        break;
                if((n-nn)>1)
                    mexErrMsgIdAndTxt("eig_omp:defectivematrix","Eigenvectors of defective eigenvalues cannot be orthogonalised.");
                // Re-interleave to Matlab output complex matrix format
                n = 0;
                for(jj=ii; jj<kk; jj++) {
                    for(nn=0; nn<m; nn++) {
                        Vr[id[jj]*m+nn] = pVz[n*2*m+nn*2];
                        Vi[id[jj]*m+nn] = pVz[n*2*m+nn*2+1];
                    }
                    n++;
                }
                ii = jj-1;
            }
    }
}

// This is called for complex general matrices - both eigenvalues and eigenvectors are actually complex
void orth(mwSignedIndex m, double *D, double *V, double *work, bool isreal)
{
    int *id, ii, jj, kk, nn;
    mwSignedIndex n, info;
    double *Dr, *pVz, *pS;
    double abstol = sqrt(DBL_EPSILON);
    char jobu = 'O';
    char jobvt = 'N';
    mwSignedIndex lwork = 3*m;
    double *zwork = work + (2*m*(m+7)) - 2*3*m;
    double *rwork = zwork - 5*m;
    size_t msz = 2*m*sizeof(double);

    // Assume workspace is size [m*(m+7)]*sizeof(complexdouble)
    id = (int*)work;
    Dr = work + m;
    for(ii=0; ii<(int)m; ii++) {
        id[ii] = ii;
        Dr[ii] = D[2*ii];
    }

    // Use QuickSort to get indices of eigenvalues in ascending order by real part
    quicksort(id, Dr, m);

    pVz = work+2*m;
    for(ii=0; ii<(int)m; ii++) {
        if((ii+1)<(int)m)
            if(fabs(Dr[id[ii]]-Dr[id[ii+1]])<abstol) {
                kk = ii;
                while(true)
                    if(fabs(Dr[id[ii]]-Dr[id[++kk]])>abstol || kk>=(int)(m-1)) 
                        break;
                // Populates a complex matrix with the degenerate eigenvectors, taking care about interleaving
                n = 0;
                for(jj=ii; jj<kk; jj++)
                    memcpy(&pVz[(n++)*2*m], &V[id[jj]*2*m], msz);
                // Does the singular value decomposition to get the orthogonal basis and singular values
                pS = work + 2*(n+1)*m;
                zgesvd(&jobu, &jobvt, &m, &n, pVz, &m, pS, pVz, &m, pVz, &m, zwork, &lwork, rwork, &info);
                // Checks that number of singular values == n
                for(nn=0; nn<n; nn++)
                    if(fabs(pS[nn])<0)
                        break;
                if((n-nn)>1)
                    mexErrMsgIdAndTxt("eig_omp:defectivematrix","Eigenvectors of defective eigenvalues cannot be orthogonalised.");
                // Permutes back
                n = 0;
                for(jj=ii; jj<kk; jj++)
                    memcpy(&V[id[jj]*2*m], &pVz[(n++)*2*m], msz);
                ii = jj;
            }
    }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    mwSignedIndex m, n, nd;
    const size_t *dims;
    int ib, nblock, nb, do_sort=0;
    int *blkid, valint;
    char jobz, *parstr, *valstr;
    bool *issym, anynonsym=false, do_orth=false;
    // Tolerance on whether matrix is symmetric/hermitian
    double tolsymm = sqrt(DBL_EPSILON);
    int nthread = omp_get_max_threads();
//  mexPrintf("Number of threads = %d\n",nthread);

    // Checks inputs
//  if(nrhs!=1) {
//      mexErrMsgIdAndTxt("eig_omp:nargin","Number of input argument must be 1.");
//  }
    if(!mxIsNumeric(prhs[0])) {
        mexErrMsgIdAndTxt("eig_omp:notnumeric","Input matrix must be a numeric array.");
    }
    nd = mxGetNumberOfDimensions(prhs[0]);
    if(nd<2 || nd>3) {
        mexErrMsgIdAndTxt("eig_omp:wrongdims","Only 2D or 3D arrays are supported.");
    }
    dims = mxGetDimensions(prhs[0]);
    m = dims[0];
    n = dims[1];
    if(m!=n) {
        mexErrMsgIdAndTxt("eig_omp:notsquare","Input matrix is not square.");
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
            if(strcmp(parstr,"orth")==0) {
                if((nrhs-1)>ib) {
                    if(!mxIsChar(prhs[ib+1]) && (mxGetM(prhs[ib+1])*mxGetN(prhs[ib+1]))>0) {
                        if(mxIsLogical(prhs[ib+1]))
                            valint = *(int*)mxGetData(prhs[ib+1]);
                        else 
                            valint = (int)(*(double*)mxGetData(prhs[ib+1]));
                        do_orth = (valint!=0) ? true : false;
                        ib++;
                    }
                    else
                        do_orth = true;
                }
                else {
                    do_orth = true;
                }
            }
            else if(strcmp(parstr,"sort")==0) {
                if((nrhs-1)>ib) {
                    if(mxIsChar(prhs[ib+1])) {
                        valstr = mxArrayToString(prhs[ib+1]);
                        if(strcmp(valstr,"orth")==0)
                            continue;
                        else if(strcmp(valstr,"ascend")==0)
                            do_sort = 1;
                        else if(strcmp(valstr,"descend")==0)
                            do_sort = -1;
                        else 
                            mexErrMsgIdAndTxt("eig_omp:badsortarg","Arguments to 'sort' keyword must be either 'ascend' or 'descend'.");
                    }
                    else {
                        if(mxIsLogical(prhs[ib+1]))
                            valint = *(int*)mxGetData(prhs[ib+1]);
                        else 
                            valint = (int)(*(double*)mxGetData(prhs[ib+1]));
                        do_sort = (valint==-1||valint==0) ? valint : 1;
                    }
                    ib++;
                }
                else
                    do_sort = 1;
            }
        }
    }
//  mexPrintf("do_orth = %d; do_sort = %d\n",do_orth,do_sort);

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

#pragma omp parallel default(none) shared(plhs,prhs) \
    firstprivate(nthread, m, nlhs, nd, ib, blkid, jobz, anynonsym, issym, do_orth, do_sort)
    {
#pragma omp for
        for(int nt=0; nt<nthread; nt++) {
            // These variables must be declared within the loop to make them local (and private)
            //   or we get memory errors.
            double *M, *E, *V=0, *D, *Di, *ptr_M, *ptr_Mi, *ptr_V, *ptr_Vi;
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
            lwork = do_orth ? ( (26*m>(2*m*(m+7))) ? 26*m : (2*m*(m+7)) )  
                            : ( (26*m>(2*m*(m+3))) ? 26*m : (2*m*(m+3)) );
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
                    if(anynonsym) {
                        zgeev(&jobz, &jobzn, &m, M, &lda, D, V, &ldz, V, &ldz, zwork, &lzwork, work, &info);
                        if(do_sort)
                            if(nlhs>1)
                                sort(m, D, V, work, do_sort);
                            else 
                                sort(m, D, 0, work, do_sort);
                        if(do_orth)
                            orth(m, D, V, work, 0);
                    }
                    else {
                        zheevr(&jobz, &range, &uplo, &m, M, &lda, &vl, &vu, &il, &iu, &abstol, &numfind, 
                                D, V, &ldz, isuppz, zwork, &lzwork, work, &lwork, iwork, &liwork, &info);
                        // ZHEEVR outputs eigenvectors in ascending order by default.
                        if(do_sort==-1) {
                            fliplr(D,1,m,work,1);
                            if(nlhs>1) 
                                fliplr(V,m,m,work,0);
                        }
                    }
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
                    if(anynonsym) {
                        dgeev(&jobzn, &jobz, &m, M, &lda, D, Di, V, &ldz, V, &ldz, work, &lwork, &info);
                        if(do_sort)
                            if(nlhs>1)
                                sort(m, D, Di, V, work, do_sort);
                            else
                                sort(m, D, Di, 0, work, do_sort);
                        if(nlhs>1 && do_orth)
                            orth(m, D, Di, V, mxGetPi(plhs[0])+ib*m2, work, 1);
                    }
                    else {
                        dsyevr(&jobz, &range, &uplo, &m, M, &lda, &vl, &vu, &il, &iu, &abstol, &numfind, 
                                D, V, &ldz, isuppz, work, &lwork, iwork, &liwork, &info);
                        // DSYEVR outputs eigenvectors in ascending order by default.
                        if(do_sort==-1) {
                            fliplr(D,1,m,work,1);
                            if(nlhs>1) 
                                fliplr(V,m,m,work,1);
                        }
                    }
                    // Need to account for complex conjugate eigenvalue/vector pairs (imaginary parts 
                    //   of eigenvectors stored in consecutive columns)
                    if(nlhs>1 && !issym[ib] && !do_orth) {
                        ptr_V = mxGetPr(plhs[0]) + ib*m2;
                        ptr_Vi = mxGetPi(plhs[0]) + ib*m2;
                        // Complex conjugate pairs of eigenvalues appear consecutively with the eigenvalue
                        //   having the positive imaginary part first.
                        for(ii=0; ii<m; ii++) {
                            if(*(Di+ii)>0.) {
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
