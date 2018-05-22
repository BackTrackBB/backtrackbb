/* Routines to diagonalize a 3x3 upper triangular matrix
 */
#include "IA_Kdiag.h"

/* use Charlier et al. approximation or not ? */
#define APX 0

/* matrix by matrix multiplication */
inline void mat_mult(double **A, double **B, int a_rows, int a_columns, int b_columns, double **C)
{
    double sum=0.0;
    int i,k,l;

    /*sanity check */
    if( (A == NULL) || (B == NULL) || (C== NULL)) {
        fprintf(stderr, "mat_mult NULL pointer for input/output matrices");
        exit(1);
    }

    for(l=0; l<b_columns; l++) {
        for(i=0; i<a_rows; i++) {
            sum = 0.0;
            for(k=0; k<a_columns; k++) {
                sum += A[i][k]*B[k][l];
            }
            C[i][l]=sum;
        }
    }
}


inline double **mat_transp(double **A,  int n_rows, int n_columns)
{
    double **B=ealloc2double(n_rows,n_columns);
    int i,j;

    for(i=0; i<n_columns; i++) {
        for(j=0; j<n_rows; j++) {
            B[i][j]=A[j][i];
        }
    }
    return(B);
}


/* compute norm squared */
#ifdef _MSC_VER
double norm2(double *a, int n_rows)
#else
inline double norm2(double *a, int n_rows)
#endif
{
    int i;
    double sum=0.0;

    for(i=0; i<n_rows; i++)
        sum+=a[i]*a[i];

    return(sum);
}


void KogExact(double a00, double a01, double a11, double **Q, double **C, double **E)
{
    /*  -----------------------------------------------------------------------
     * Exact variant of Kogbetlianz's algorithm after
     *
     * Charlier, J. P. and Vanbegin, M. and Dooren, P. Van, "On Efficient
     * Implementations of Kogbetliantz's Algorithm for Computing the
     * Singular Value Decomposition", Numerische Mathematik,
     * vol. 52, pp. 279-300, 1988
     *
     * Computes left and right rotation matrices Q and E as well as
     * diagonal matrix A from upper triangular 2x2 matrix with elements
     * a_ij
     * -----------------------------------------------------------------------*/


    /*sanity check */
    if( (Q ==NULL) || (C == NULL) || (E==NULL)) {
        fprintf(stderr, "Null pointer for output matrices");
        exit(1);
    }

    double sigma;
    double t_phi=0.0;
    double t_psi=0.0;
    double C_phi,S_phi;
    double C_psi,S_psi;
    double f,g,x;

    int u1,u2;

    double a = a11;
    double b = a00;
    double c = a01;

    memset(*Q,0,4*sizeof(double));
    Q[0][0]=Q[1][1]=1.0;
    memset(*C,0,4*sizeof(double));
    memset(*E,0,4*sizeof(double));
    E[0][0]=E[1][1]=1.0;

    /* Return -almost- untouched if extra diagonal element already small ~ 0 */

    C[0][0]=b;
    C[1][1]=a;
    C[0][1]=0.0;
    C[1][0]=0.0;


    if ( fabs(c) > KOG_EPS ) {

        if( fabs(a) > fabs(b)) {

            sigma = 2*b*c/((a-b)*(a+b)+c*c); // if(b == 0) sigma=0
            t_phi = sigma/(1+sqrt(1+sigma*sigma)); // if ( b==0 ) t_phi = 0

            if ( a != 0.0) {
                t_psi = (c+t_phi*b)/a;
            }

        }
        else {

            sigma = 2*a*c/((b-a)*(b+a)+c*c);  // if(a == 0) sigma=0
            t_psi = -sigma/(1+sqrt(1+sigma*sigma)); // if(a == 0) t_psi=0

            if( b != 0.0) {
                t_phi = (t_psi*a-c)/b;
            }
        }

        C_phi = 1.0/sqrt(1.0+t_phi*t_phi);
        S_phi = t_phi*C_phi;
        C_psi = 1.0/sqrt(1+t_psi*t_psi);
        S_psi = t_psi*C_psi;
        x = C_psi/C_phi;

        f = b*x;
        C[0][0] = fabs(f);

        u1 = (f < 0.0) ? -1 : 1;

        g = a/x;
        C[1][1] = fabs(g);

        u2 = (g < 0.0) ? -1 : 1;

        if( fabs(g) <= fabs(f)) {
            Q[0][0] = u1*C_psi;
            Q[0][1] = u2*S_psi;
            Q[1][0] = -u1*S_psi;
            Q[1][1] = u2*C_psi;

            E[0][0] = C_phi;
            E[0][1] = S_phi;
            E[1][0] =-S_phi;
            E[1][1] = C_phi;
        }
        else {

            /* sort singular values descending : permute Q,A and E
             */

            /* switch columns on Q */
            Q[0][0] = u2*S_psi;
            Q[0][1] = u1*C_psi;
            Q[1][0] = u2*C_psi;
            Q[1][1] = -u1*S_psi;


            /* switch coumns on E */

            E[0][0] = S_phi;
            E[0][1] = C_phi;
            E[1][0] = C_phi;
            E[1][1] =-S_phi;

            /* permute A : switch rows and columns */

            f = C[0][0];
            C[0][0] = C[1][1];
            C[1][1] = f;

        }

    } /*  if ( fabs(c) > 0.0 ) */
}


void KogApx(double a00, double a01, double a11, double **Q, double **C, double **E)
{
    /*  -----------------------------------------------------------------------
     * Approximation variant of Kogbetlianz's algorithm after
     *
     * Charlier, J. P. and Vanbegin, M. and Dooren, P. Van, "On Efficient
     * Implementations of Kogbetliantz's Algorithm for Computing the
     * Singular Value Decomposition", Numerische Mathematik,
     * vol. 52, pp. 279-300, 1988
     *
     * Computes left and right rotation matrices Q and E as well as
     * diagonal matrix A from upper triangular 2x2 matrix with elements
     * a_ij
     * -----------------------------------------------------------------------*/

    /*sanity check */
    if( (Q ==NULL) || (C == NULL) || (E==NULL)) {
        fprintf(stderr, "Null pointer for output matrices");
        exit(1);
    }

    double sigma;
    double t_phi=0.0;
    double t_psi=0.0;
    double C_phi,S_phi;
    double C_psi,S_psi;
    double f,g,x;

    int u1,u2;

    double a = a11;
    double b = a00;
    double c = a01;

    memset(*Q,0,4*sizeof(double));
    Q[0][0]=Q[1][1]=1.0;
    memset(*C,0,4*sizeof(double));
    memset(*E,0,4*sizeof(double));
    E[0][0]=E[1][1]=1.0;

    /* Return -almost- untouched if diagonal element already small ~ 0 */

    C[0][0]=b;
    C[1][1]=a;
    C[0][1]=0.0;
    C[1][0]=0.0;

    if ( fabs(c) > KOG_EPS ) { /* there is interaction with termination condition! */

        if( fabs(a) > fabs(b)) {

            sigma = b*c/((a-b)*(a+b)+c*c); // if(b == 0) sigma=0
            t_phi = sigma/(1+sigma*sigma);

            if( a != 0.0 ) {
                t_psi = (c+t_phi*b)/a;
            }
        }
        else {

            sigma = a*c/((b-a)*(b+a)+c*c); // if(a == 0) sigma=0
            t_psi = -sigma/(1+sigma*sigma);

            if( b != 0.0 ) {
                t_phi = (t_psi*a-c)/b;
            }
        }

        C_phi = 1.0/sqrt(1.0+t_phi*t_phi);
        S_phi = t_phi*C_phi;
        C_psi = 1.0/sqrt(1+t_psi*t_psi);
        S_psi = t_psi*C_psi;
        x = C_psi/C_phi;

        f = b*x;
        C[0][0] = fabs(f);


        u1 = (f < 0.0) ? -1 : 1;

        g = a/x;
        C[1][1] = fabs(g);

        u2 = (g < 0.0) ? -1 : 1;

        if( fabs(g) <= fabs(f)) {
            Q[0][0] = u1*C_psi;
            Q[0][1] = u2*S_psi;
            Q[1][0] = -u1*S_psi;
            Q[1][1] = u2*C_psi;

            E[0][0] = C_phi;
            E[0][1] = S_phi;
            E[1][0] =-S_phi;
            E[1][1] = C_phi;
        }
        else {

            /* sort singular values descending : permute Q,A and E
             */

            /* switch columns on Q */
            Q[0][0] = u2*S_psi;
            Q[0][1] = u1*C_psi;
            Q[1][0] = u2*C_psi;
            Q[1][1] = -u1*S_psi;


            /* switch columns on E */

            E[0][0] = S_phi;
            E[0][1] = C_phi;
            E[1][0] = C_phi;
            E[1][1] =-S_phi;

            /* permute C : switch rows and columns */

            f = C[0][0];
            C[0][0] = C[1][1];
            C[1][1] = f;

        }

    } /*  if ( fabs(c) > 0.0 ) */
}


void KDiag2by2(double **A, int i, int j, double **R, double **Y, double **T, char apx)
{
    /* -----------------------------------------------------------------------
       Diagonalize 2x2 upper triangular sub-matrix A_ii,A_ij,A_jj of 3x3 matrix
       A; return A's submatrix SVD in 3x3 matrices R,Y,T.
       Use approximate Kb algorithm if apx == 1
       R,Y and T must be pre-allocated
       -----------------------------------------------------------------------*/

    /* sanity */
    if (( A == NULL ) || (R==NULL) || (Y==NULL) || (T==NULL)) {
        fprintf(stderr, "Nullpointer in KDiag2by2");
        exit(1);
    }
    if(j <= i) {
      fprintf(stderr, "j > i !");
      exit(1);
    }

    double **Q = ealloc2double(2,2);
    double **C = ealloc2double(2,2);
    double **E = ealloc2double(2,2);
    double **TMP= ealloc2double(3,3);


    if ( apx ) {
        KogApx(A[i][i], A[i][j], A[j][j], Q, C, E);
    }
    else {
        KogExact(A[i][i], A[i][j], A[j][j], Q, C, E);
    }

    R[0][0]=R[1][1]=R[2][2]=1.0;
    R[0][1]=R[0][2]=R[1][0]=R[1][2]=R[2][0]=R[2][1]=0.0;

    T[0][0]=T[1][1]=T[2][2]=1.0;
    T[0][1]=T[0][2]=T[1][0]=T[1][2]=T[2][0]=T[2][1]=0.0;

    R[i][i]=Q[0][0];
    R[i][j]=Q[0][1];
    R[j][i]=Q[1][0];
    R[j][j]=Q[1][1];

    T[i][i]=E[0][0];
    T[i][j]=E[0][1];
    T[j][i]=E[1][0];
    T[j][j]=E[1][1];

    /* rotate A; Y=R'*A*T */

    double **RT = mat_transp(R,3,3);

    mat_mult(RT, A, 3, 3, 3,TMP);  /* optimize here for speed */
    mat_mult(TMP, T, 3, 3, 3,Y);

    /* prevent error accumulation on lower triangle, force to zero */

    Y[1][0]=0.0;
    Y[2][0]=0.0;
    Y[2][1]=0.0;

    free2double(RT);
    free2double(TMP);
    free2double(Q);
    free2double(C);
    free2double(E);
}


void BA3_SVD(double **A, double **U, double **S, double **V, int *n)
{
    /*  -----------------------------------------------------------------------
     * Singular value Decomposition of a 3 x 3 upper triangular matrix A.
     * On output n is number of iterations
     * Matrix A is destroyed !!
     *  ----------------------------------------------------------------------- */

    double keps;
    char apx;

    double **R=ealloc2double(3,3);
    double **Y=ealloc2double(3,3);
    double **T=ealloc2double(3,3);
    double **TMP=ealloc2double(3,3);

    double p;
    int i,j,k;

    if( (U == NULL) || (Y == NULL) || (V ==NULL)) {
        fprintf(stderr, "Null pointer(s) in BA3_SVD");
        exit(1);
    }

    memset(*U,0,9*sizeof(double));
    memset(*V,0,9*sizeof(double));

    U[0][0]=U[1][1]=U[2][2]=1.0;
    V[0][0]=V[1][1]=V[2][2]=1.0;

    *n=0;

    /* iterate upper sub-diagonal */

    apx=APX;

    do {

        KDiag2by2(A,1,2,R,Y,T,apx);  /* 1,2 element A -> Y*/

        mat_mult(U,R,3,3,3,TMP);
        memcpy(*U,*TMP,9*sizeof(double)); // free(TMP)??

        mat_mult(V,T,3,3,3,TMP);
        memcpy(*V,*TMP,9*sizeof(double));


        KDiag2by2(Y,0,1,R,A,T,apx); /* 0,1 element  Y - >A */

        mat_mult(U,R,3,3,3,TMP);
        memcpy(*U,*TMP,9*sizeof(double));

        mat_mult(V,T,3,3,3,TMP);
        memcpy(*V,*TMP,9*sizeof(double));


        (*n)++;
        keps = A[0][1]*A[0][1] + A[1][2]*A[1][2]; /* norm^2 of off-diagonal */

        if(keps < KOG_EPS) apx=0;

    } while( (keps > (KOG_EPS * KOG_EPS) ) && (*n < KOG_MAXIT) );



    /* eliminate upper corner element 0,2 */

    KDiag2by2(A,0,2,R,Y,T,0); /* A - > Y */

    mat_mult(U,R,3,3,3,TMP);
    memcpy(*U,*TMP,9*sizeof(double)); // U=U*R

    mat_mult(V,T,3,3,3,TMP);
    memcpy(*V,*TMP,9*sizeof(double)); // V=V*T

    memset(*S,0,9*sizeof(double));

    S[0][0]=Y[0][0];
    S[1][1]=Y[1][1];
    S[2][2]=Y[2][2];

    /*FIXED: do permutation so that S0 > S1 >S2 !!! */

    for (i=0; i<3; i++) {
        k=i;
        p=S[k][k];
        for (j=i+1; j<3; j++)
            if (S[j][j] > p) {
                k=j;
                p=S[k][k];
            }
        if (k != i) {
            S[k][k]=S[i][i];
            S[i][i]=p;
            for (j=0; j<3; j++) {
                p=V[i][j];
                V[i][j]=V[k][j];
                V[k][j]=p;
            }
            for (j=0; j<3; j++) {
                p=U[j][i];
                U[j][i]=U[j][k];
                U[j][k]=p;
            }
        }
    }

    free2double(R);
    free2double(Y);
    free2double(T);
    free2double(TMP);
}


#ifdef IA_DEBUG

void printvect(double *d, int m)
{
    int i;

    fprintf(stderr,"\n\n");

    for(i=0; i<m; i++) {
        fprintf(stderr,"%f\n",d[i]);
    }
    fprintf(stderr,"\n\n");
}


void printmat(double **M, int m, int n)
{
    int i,j;

    fprintf(stderr,"\n\n");
    for(i=0; i<m; i++) {
        for (j=0; j<n; j++) {
            fprintf(stderr,"%e ",M[i][j]);
        }
        fprintf(stderr,"\n");
    }
    fprintf(stderr,"\n\n");
}


void printUSV(double **U, double *S, double **V, int m, int n)
{
    int i,j,k;
    double **tmp;

    tmp= ealloc2double(n,m);

    /* tmp=U*S */
    for(i=0; i<m; i++) {
        for(k=0; k<n; k++) {
            tmp[i][k]+=U[i][k]*S[k];

        }
    }

    printmat(tmp,m,n);

    /* U= tmp * V' */

    for(i=0; i<m; i++) {
        for(k=0; k<n; k++) {
            U[i][k]=0;
            for(j=0; j<n; j++) {
                U[i][k]+=tmp[i][j]*V[k][j];
            }
        }
    }
    printmat(U,m,n);

    free2double(tmp);

}
#endif /*IA_DEBUG */
