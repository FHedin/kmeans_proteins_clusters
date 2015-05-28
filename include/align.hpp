/*
 * kmeans : c++ implementation of the kmeans algorithm (http://en.wikipedia.org/wiki/K-means_clustering)
 * for finding clusters of stable microstates when studying ligand migration in proteins
 * 
 * Copyright (c) 2015, Florent Hédin, Pierre-André Cazade, and the University of Basel.
 * All rights reserved.
 * 
 * The 3-clause BSD license is applied to this software.
 * See LICENSE.txt
 */

#ifndef ALIGN_HPP_INCLUDED
#define ALIGN_HPP_INCLUDED
    
namespace align
{

typedef struct
{
    double X;
    double Y;
    double Z;
} POINT_3D;

#define LAPACK_SGESVD_ERROR             120
#define SIGN(a)                (a==0.0)?0.0:(a/fabs(a))
    
void kabsch_align(std::vector<float>& x, std::vector<float>& y, std::vector<float>& z,
             std::vector<float>& xr, std::vector<float>& yr, std::vector<float>& zr,
             std::vector<bool>& index_c);

POINT_3D getCentroid(std::vector<float>& x, std::vector<float>& y, std::vector<float>& z, std::vector<bool>& index_c);

void center_sys(std::vector<float>& x, std::vector<float>& y, std::vector<float>& z, std::vector<bool>& index_c);

void get_cov_mat(std::vector<float>& x, std::vector<float>& y, std::vector<float>& z,
                 std::vector<float>& xr, std::vector<float>& yr, std::vector<float>& zr,
                 std::vector<bool>& index_c, float covmat[9]);

double determinant_3x3(float inp[9]);

void matmul_3x3(float a[9], float b[9], float c[9]);

void svd_f(float a[9], float u[9], float vt[9]);

void rot_mat(double d, float u[9], float vt[9], std::vector<float>& x, std::vector<float>& y, std::vector<float>& z);

#ifdef  __cplusplus
extern "C"
{
#endif
    
// for the lapack routine sgesvd : singular values decomposition fo floats
extern void sgesvd_(char *jobu, char *jobvt, int *m, int *n, float *a, int *lda, float *s, float *u, int *ldu,
                    float *vt, int *ldvt, float *work, int *lwork, int *info);
//for the blas routine dgemm : matrix multiplication for floats
extern void sgemm_(char* transA, char* transB, int *M, int *N, int * K, float *alpha, float *A, int *ldA, float *B, int *ldB, float *beta, float *C, int *ldC);

#ifdef  __cplusplus
}
#endif

}

#endif // ALIGN_HPP_INCLUDED
