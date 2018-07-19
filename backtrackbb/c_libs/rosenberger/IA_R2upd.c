/*
 * rank 2 update SVD, see Rosenberger, BSSA, June 2010
 *
 * Revision 1.5  2018/05/22 Claudio Satriano <satriano@ipgp.fr>
 *               - compiling under Win10
 *
 * Revision 1.4  2016/03/15 Claudio Satriano <satriano@ipgp.fr>
 *               - subtrack2 now only returns the polarization filter
 *
 * Revision 1.3  2014/05/11 Claudio Satriano <satriano@ipgp.fr>
 *               - Added rl_filter option to subtrack2
 *
 * Revision 1.2  2011/01/31 22:59:52  andreas
 * Added param char proj to subtrack2.
 * 0<=proj<=2 for KLT projection using 1 or 2 eigen-images
 * proj = 0 turns KLT off, original data are returned
 *
 * Revision 1.1  2011/01/31 22:57:05  andreas
 * Initial revision
 *
 */
#include "IA_Kdiag.h"
//#define DEBUG


USV32_struct *New_USV32_struct()
{
    USV32_struct *P;

    if((P=(USV32_struct *)malloc(sizeof(USV32_struct))) == NULL) {
        fprintf(stderr, "Out of mem in New_USV32_struct()\n");
        exit(1);
    }

    P->U=ealloc2double(2,3);
    P->S=ealloc2double(2,2);

    return(P);
}


void Destroy_USV32_struct(USV32_struct *P)
{
    if(P!=NULL) {
        if(P->U!=NULL) free2double(P->U);
        if(P->S!=NULL) free2double(P->S);
        free(P);
    }
}


void subtrack2(const IA_Dvect *x, double *pol_out,
               const float lambda, const float delta,
               char proj, char rl_filter,
               USV32_struct *P, char r2u_init)
{
    /* Input:  3-C data-vector x
     * Output: 3-C data-vector with P- type polarization: p
     * Output: 3-C data-vector with S- type polarization: s
     *
     * lambda: exponential forgetting factor 0<lambda<=1.0
     *
     * delta:  minumum innovation; machine_EPS << delta <= instrument sensitivity
     *
     * P:   3x2 Workspace matrices U, S, calibration values for Z, N, E channels
     */
    double **Ut = ealloc2double(3,3);
    double **St = ealloc2double(3,3);
    double **Vt = ealloc2double(3,3);
    double **K = ealloc2double(3,3);
    double **Temp = ealloc2double(2,3);

    double m[2];
    double p[3];
    double y[3];
    double z[3];

    double r;
    double norm_p;
    double r_ortho;

    double pol; /* incident angle = acos(U[0][0]) */
    double rl;  /* rectilinearity = 1 - S[1][1]/S[0][0] */

    int i,j,k,n=0;


    /* sanity */
    if (P == NULL) {
        fprintf(stderr, "Subtrack2: Nullvector P\n");
        exit(1);
    }

    if((proj < 0) || (proj > 2)) {
        fprintf(stderr, "Subtrack2: 0 <= proj <= 2\n");
        exit(1);
    }

    /* convert input  to double  */
    z[0] = x->z;
    z[1] = x->n;
    z[2] = x->e;

    /* Initialization ------------------------------------------ */
    if (r2u_init) {
        P->U[0][0] = 1.0;
        r = sqrt(norm2(z,3));
        if (r > 0) {
            for (i=0; i<3; i++) {
                P->U[i][0] = z[i]/r;
                P->U[i][1] = 0.0;
            }
        }
        P->S[0][0] = r;
    }
    /* End initialization ---------------------------------------
     */

    /* project on current subspace: m = U'x */
    for (i=0; i<2; i++) {
        m[i] = 0.0;
        for (j=0; j<3; j++) {
            m[i] += P->U[j][i]*z[j];
        }
    }

    /* auxiliary vector: y = Um */
    for (i=0; i<3; i++) {
        y[i] = 0.0;
        for (j=0; j<2; j++) {
            y[i] += P->U[i][j]*m[j];
        }
    }

    /* project on orthogonal subspace p = x - Um */
    for (i=0; i<3; i++) {
        p[i] = z[i]-y[i];
    }

    /*
     * |p|=sqrt(x'x -2m'm + (Um)'Um), for better fp accuracy
     */
    r = norm2(z, 3) - 2*norm2(m, 2) + norm2(y, 3);

    if (r < 0.0) { /* must be small */
        if (r < -delta) fprintf(stderr, "subtrack2: negative r= %e\n", r);
        norm_p = 0.0;
    } else {
        norm_p = sqrt(r);
    }

    /* rank increase or not ? */
    norm_p = norm_p > delta ? norm_p : 0.0;

    /* construct update matrix K */
    memset(*K, 0, 9*sizeof(double));
    K[0][0] = (1-lambda)*P->S[0][0];
    K[1][1] = (1-lambda)*P->S[1][1];
    K[0][2] = m[0];
    K[1][2] = m[1];
    K[2][2] = norm_p;

    /* Compute SVD of K */
    BA3_SVD(K, Ut, St, Vt, &n); /*  Kogbetlianz' SVD */
    //printf("%f %f %f\n", St[0][0], St[1][1], St[2][2]);

#ifdef DEBUG
    if (St[0][0]*St[1][1]*St[2][2] <= 0)
        fprintf(stderr,"Neg. or Zero  Sing Val \n");
#endif

    /* discard smallest singular value */
    P->S[0][0] = St[0][0];
    P->S[1][1] = St[1][1];

    if (norm_p > 0.0) {
        /* rank reduction
         * U=[U,m_p/p]*Us(1:3,1:2);
         */
        for (i=0; i<3; i++) {
            for (j=0; j<2; j++) {
                Temp[i][j] = 0.0;
                for (k=0; k<2; k++) {
                    Temp[i][j] += P->U[i][k]*Ut[k][j];
                }
                /* contribution from virtual 3rd column of P->U */
                Temp[i][j] += p[i]/norm_p*Ut[2][j];
            }
        }
    }
    else {
        /* rank preserving P->U 3 x 2, Ut 3 x 3
         * U=U*Ut(1:2,1:2)
         */
        for (i=0; i<3; i++) {
            for (j=0; j<2; j++) {
                Temp[i][j] = 0.0;
                for (k=0; k<2; k++) {
                    Temp[i][j] += P->U[i][k]*Ut[k][j];
                }
            }
        }
    }
    memcpy(*(P->U), *Temp, 6*sizeof(double));

    /* perform Ortho-Check here */
    r_ortho = 0.0;
    for (i=0; i<3; i++) {
        r_ortho += P->U[i][0]*P->U[i][1];
    }
    if (fabs(r_ortho) > 1.0e-5) {
        fprintf(stderr, "subtrack2: Orthonormality lost %e\n", r_ortho);
        /* Do something about it */
    }

    /* output */
    if (proj) {
        /* project: x = U*U'*z = U*m */
        for (i=0; i<proj; i++) {
            m[i] = 0.0;
            for (j=0; j<3; j++) {
                m[i] += P->U[j][i]*z[j];
            }
        }
        /* project on current subspace: m = U'x */
        for (i=0; i<3; i++) {
            z[i] = 0.0;
            for (j=0; j<proj; j++) {
                z[i] += P->U[i][j]*m[j];
            }
        }
    }

    if (rl_filter) {
        /* linear polarization (rectilinearity) as ratio
        *  of second over first singular value
        */
        if(P->S[0][0] < 1e-16) {
            rl = 1.0;
            fprintf(stderr,
                    "subtrack2: (Near) zero principal singular value, input=%e,%e,%e, ratio=%e\n",
                    z[0], z[1], z[2], P->S[0][0]);
        }
        else {
            rl = 1.0 - P->S[1][1]/P->S[0][0];
        }
    } else {
        rl = 1.;
    }

    pol = fabs(P->U[0][0]); /* acos(U[0][0]) is vertical angle of incidence */
    *pol_out = pol;


    free2double(Ut);
    free2double(St);
    free2double(Vt);
    free2double(K);
    free2double(Temp);
}
