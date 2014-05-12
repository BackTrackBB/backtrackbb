/* rosenberger.c
 *
 * Glue code to call subtrack2() on a three-component time series
 *
 * (c) 2014 - Claudio Satriano <satriano@ipgp.fr>
 *
 */
#include "IA_Kdiag.h"

void rosenberger(double *dataX, double *dataY, double *dataZ,
                 double *dataX_P, double *dataY_P, double *dataZ_P,
                 double *dataX_S, double *dataY_S, double *dataZ_S,
                 int npts,
                 float lambda, float delta, char proj, char rl_filter)
{
    int n;
    IA_Dvect data;
    IA_Dvect dataP;
    IA_Dvect dataS;
    USV32_struct *P = New_USV32_struct();
    char r2u_init = 0x1;

    for (n=0; n<npts; n++) {
        data.e = dataX[n];
        data.n = dataY[n];
        data.z = dataZ[n];

        subtrack2(&data, &dataP, &dataS, lambda, delta,
                  proj, rl_filter, P, r2u_init);
        r2u_init = 0x0;

        dataX_P[n] = dataP.e;
        dataY_P[n] = dataP.n;
        dataZ_P[n] = dataP.z;
        dataX_S[n] = dataS.e;
        dataY_S[n] = dataS.n;
        dataZ_S[n] = dataS.z;
    }

    Destroy_USV32_struct(P);
}
