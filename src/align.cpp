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

#include <vector>

#include <cstdlib>
#include <cstdio>
#include <cmath>

#include "align.hpp"

using namespace std;


/**
 * Aligning a set of coordinates to a reference one using the Kabsch Algorithm :
 * http://en.wikipedia.org/wiki/Kabsch_algorithm
 */
namespace ALIGN {

// align system in x y z using references from xr yr zr
void align_to_ref(vector<float>& x, vector<float>& y, vector<float>& z,
             vector<float>& xr, vector<float>& yr, vector<float>& zr,
             vector<bool>& index_c)
{
    //covariance matrix 3x3
    float covmat[3 * 3] = { 0.f };
    //matrices 3x3 for SVD
    float u[3 * 3] = { 0.f }, vt[3 * 3] = { 0.f };

    center_sys(x, y, z, index_c);

    get_cov_mat(x, y, z,
                xr, yr, zr,
                index_c, covmat);

    double d = determinant_3x3(covmat);
    d = SIGN(d);

    svd_f(covmat, u, vt);

    rot_mat(d, u, vt, x, y, z);
}

//get the centroid of a system
POINT_3D getCentroid(vector<float>& x, vector<float>& y, vector<float>& z, vector<bool>& index_c)
{
    uint i = 0, k = 0;
    POINT_3D c = { 0., 0., 0. };

    for ( i = 0; i < x.size(); i++ )
    {
        if ( index_c[i] )
        {
            c.X += x[i];
            c.Y += y[i];
            c.Z += z[i];
            k++;
        }
    }

    c.X /= ( double ) k;
    c.Y /= ( double ) k;
    c.Z /= ( double ) k;

    return c;
}

//center the system so that the centroid is at the origin
void center_sys(vector<float>& x, vector<float>& y, vector<float>& z, vector<bool>& index_c)
{
    uint i = 0;
    POINT_3D c = { 0., 0., 0. };

    c = getCentroid(x, y, z, index_c);

    for ( i = 0; i < x.size(); i++ )
    {
        x[i] -= c.X;
        y[i] -= c.Y;
        z[i] -= c.Z;
    }
}

//get the covariance matrix
void get_cov_mat(vector<float>& x, vector<float>& y, vector<float>& z,
                 vector<float>& xr, vector<float>& yr, vector<float>& zr,
                 vector<bool>& index_c, float covmat[9])
{
//     int i = 0, j = 0, k = 0, l = 0;
//     for ( i = 0; i < 3; i++ )
//     {
//         for ( j = 0; j < 3; j++ )
//         {
//             covmat[i + 3 * j] = 0.f;
//             l = 0;
//             for ( k = 0; k < x.size(); k++ )
//             {
//                 if ( index_c[k] )
//                 {
//                     covmat[i + 3 * j] += crd[l][i] * ref[l][j];
//                     l++;
//                 }
//             }
//         }
//     }

    float *ptc = nullptr;
    float *ptr = nullptr;

    for (int i = 0; i < 3; i++ )
    {
        if(i==0)
            ptc = x.data();
        else if(i==1)
            ptc = y.data();
        else if(i==2)
            ptc = z.data();
        
        for (int j = 0; j < 3; j++ )
        {

            if(j==0)
                ptr = xr.data();
            else if(j==1)
                ptr = yr.data();
            else if(j==2)
                ptr = zr.data();
            
            covmat[i + 3 * j]=0.f;
            for (uint k = 0; k < x.size(); k++ )
            {
                if(index_c[k])
                    covmat[i + 3 * j] += ptc[k]*ptr[k];
            }

        }
    }

}

//get determinant of a 3x3 matrix stored in col order in an array 1d
double determinant_3x3(float inp[9])
{
    int i, j;
    double det = 0.0;
    double m[3][3] = {
        {0.0 }
    };

    for ( i = 0; i < 3; i++ )
        for ( j = 0; j < 3; j++ )
            m[i][j] = inp[i + 3 * j];

    det = m[0][0]*(m[1][1] * m[2][2] - m[1][2] * m[2][1]);
    det += -m[0][1]*(m[1][0] * m[2][2] - m[1][2] * m[2][0]);
    det += m[0][2]*(m[1][0] * m[2][1] - m[1][1] * m[2][0]);

    return det;
}

void matmul_3x3(float a[9], float b[9], float c[9])
{
    char trans_a = 'N';
    char trans_b = 'N';
    float one = 1.0;
    int dim = 3;

    sgemm_(&trans_a, &trans_b, &dim, &dim, &dim, &one, a, &dim, b, &dim, &one, c, &dim);
}

//singular value decomposition of the covariance matrix a
void svd_f(float a[9], float u[9], float vt[9])
{
    /**
     *        subroutine SGESVD   (   CHARACTER   JOBU,
     *                                CHARACTER   JOBVT,
     *                                INTEGER     M,
     *                                INTEGER     N,
     *                                REAL, dimension( lda, * )   A,
     *                                INTEGER     LDA,
     *                                REAL, dimension( * )    S,
     *                                REAL, dimension( ldu, * )   U,
     *                                INTEGER     LDU,
     *                                REAL, dimension( ldvt, * )  VT,
     *                                INTEGER     LDVT,
     *                                REAL, dimension( * )    WORK,
     *                                INTEGER     LWORK,
     *                                INTEGER     INFO
     *        )
     *
     *        extern void sgesvd_(char *jobu, char *jobvt, int *m, int *n, float *a, int *lda, float *s, float *u, int *ldu,
     *                            float *vt, int *ldvt, float *work, int *lwork, int *info);
     *
     */

    char jobu = 'A', jobvt = 'A';
    int m = 3, n = 3; //dimensions of covmat a
    int lda = 3, ldu = 3, ldvt = 3; //leading dimension : here the same m for lda, etc...
    float s[3] = { 0.f };
    float wk_siz = 0.f;
    int lwork = -1;
    int info = 0;

    //first we get the optimal size of the working array in wk_siz
    sgesvd_(&jobu, &jobvt, &m, &n, a, &lda, s, u, &ldu,
            vt, &ldvt, &wk_siz, &lwork, &info);

    if ( info != 0 )
    {
        printf("LAPACK error for SGESVD : info = %d\nCheck file %s line %d\n", info, __FILE__, __LINE__);
        exit(LAPACK_SGESVD_ERROR);
    }

    lwork = ( int ) wk_siz;
    float *work = NULL;
    work = ( float* ) malloc(lwork * sizeof *work);

    //now we launch the routine ...
    sgesvd_(&jobu, &jobvt, &m, &n, a, &lda, s, u, &ldu,
            vt, &ldvt, work, &lwork, &info);

    if ( info != 0 )
    {
        printf("LAPACK error for SGESVD : info = %d\nCheck file %s line %d\n", info, __FILE__, __LINE__);
        exit(LAPACK_SGESVD_ERROR);
    }

    free(work);
}

//compute and apply the rotation Matrix
void rot_mat(double d, float u[9], float vt[9], vector<float>& x, vector<float>& y, vector<float>& z)
{
    float identitySign_matrix[3 * 3] = {
        1.f, 0.f, 0.f,
        0.f, 1.f, 0.f,
        0.f, 0.f, ( float ) d
    };

    float v[3 * 3] = { 0.f };
    float ut[3 * 3] = { 0.f };

    float tmp1[3 * 3] = { 0.f };

    float rot[3 * 3] = { 0.f };

    uint i, j, n;

    for ( i = 0; i < 3; i++ )
    {
        for ( j = 0; j < 3; j++ )
        {
            v[i + 3 * j] = vt[j + 3 * i];
            ut[i + 3 * j] = u[j + 3 * i];
        }
    }

    matmul_3x3(v, identitySign_matrix, tmp1);
    matmul_3x3(tmp1, ut, rot);

    float *ptc = nullptr;
    for ( n = 0; n < x.size(); n++ )
    {
        float work[3] = { 0.f };

        for ( i = 0; i < 3; i++ )
        {
            for ( j = 0; j < 3; j++ )
            {
                if(j==0)
                    ptc = x.data();
                else if(j==1)
                    ptc = y.data();
                else if(j==2)
                    ptc = z.data();
                
                work[i] += rot[i + 3 * j] * ptc[n];
            }
        }

        for ( i = 0; i < 3; i++ )
        {
            if(i==0)
                ptc = x.data();
            else if(i==1)
                ptc = y.data();
            else if(i==2)
                ptc = z.data();
            
            ptc[n] = work[i];
        }
    }
}

}//end of namespace
