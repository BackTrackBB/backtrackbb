/* rosenberger.c
 *
 * Glue code to call subtrack2() on a three-component time series
 *
 * (c) 2014 - 2018 Claudio Satriano <satriano@ipgp.fr>
 *
 */
#include "IA_Kdiag.h"

void initlib_rosenberger() {}
void PyInit_lib_rosenberger() {}

#ifdef _MSC_VER
__declspec(dllexport)
#endif
void rosenberger(double *dataX, double *dataY, double *dataZ,
                 double *polarization, int npts,
                 float lambda, float delta, char proj, char rl_filter)
{
    int n;
    IA_Dvect data;
    USV32_struct *P = New_USV32_struct();
    char r2u_init = 0x1;

    for (n=0; n<npts; n++) {
        data.e = dataX[n];
        data.n = dataY[n];
        data.z = dataZ[n];

        subtrack2(&data, &polarization[n],
                  lambda, delta, proj, rl_filter, P, r2u_init);
        r2u_init = 0x0;
    }

    Destroy_USV32_struct(P);
}
