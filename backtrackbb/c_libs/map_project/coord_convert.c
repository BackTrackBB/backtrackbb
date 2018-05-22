#include "coord_convert.h"

/** function to read map transformation parameters from input line ***/
#ifdef _MSC_VER
__declspec(dllexport)
#endif
int get_transform(int n_proj, char* in_line) {
    int istat, ierr;
    double angle;

    // SDC
    double dlt1, dlt2, del, r, bc;

    /* Init constants */
    cPI = 4. * atan(1.); /* PI */
    cRPD = cPI / 180.; /* radians per degree */
    c111 = 10000.0 / 90.0; /* kilometers per degree */

    map_itype[n_proj] = MAP_TRANS_UNDEF;
    GeometryMode = MODE_RECT;

    /* read transform input line */

    sscanf(in_line, "%s", map_trans_type[n_proj]);

    if (strcmp(map_trans_type[n_proj], "GLOBAL") == 0) {

        // mode
        GeometryMode = MODE_GLOBAL;

        map_itype[n_proj] = MAP_TRANS_GLOBAL;
        istat = sscanf(in_line, "%s", map_trans_type[n_proj]);

        angle = 0.0;
        map_cosang[n_proj] = cos(angle);
        map_sinang[n_proj] = sin(angle);

        sprintf(MapProjStr[n_proj],
                "TRANSFORM  %s",
                map_trans_type[n_proj]);
        nll_putmsg(3, MapProjStr[n_proj]);

        ierr = 0;
        if (ierr < 0 || istat != 1) {
            nll_puterr("ERROR: reading GLOBAL transformation parameters");
            return (-1);
        }


    } else if (strcmp(map_trans_type[n_proj], "NONE") == 0) {

        map_itype[n_proj] = MAP_TRANS_NONE;
        istat = sscanf(in_line, "%s", map_trans_type[n_proj]);

        angle = 0.0;
        map_cosang[n_proj] = cos(angle);
        map_sinang[n_proj] = sin(angle);

        sprintf(MapProjStr[n_proj],
                "TRANSFORM  %s",
                map_trans_type[n_proj]);
        nll_putmsg(3, MapProjStr[n_proj]);

        ierr = 0;
        if (ierr < 0 || istat != 1) {
            nll_puterr("ERROR: reading NONE transformation parameters");
            return (-1);
        }


    } else if (strcmp(map_trans_type[n_proj], "SIMPLE") == 0) {

        map_itype[n_proj] = MAP_TRANS_SIMPLE;
        istat = sscanf(in_line, "%s %lf %lf %lf",
                map_trans_type[n_proj], &map_orig_lat[n_proj], &map_orig_long[n_proj],
                &map_rot[n_proj]);

        angle = -cRPD * map_rot[n_proj];
        map_cosang[n_proj] = cos(angle);
        map_sinang[n_proj] = sin(angle);

        sprintf(MapProjStr[n_proj],
                "TRANSFORM  %s LatOrig %lf  LongOrig %lf  RotCW %lf",
                map_trans_type[n_proj], map_orig_lat[n_proj], map_orig_long[n_proj], map_rot[n_proj]);
        nll_putmsg(3, MapProjStr[n_proj]);

        ierr = 0;
        if (checkRangeDouble("TRANS",
                "LatOrig", map_orig_lat[n_proj], 1, -90.0, 1, 90.0) != 0)
            ierr = -1;
        if (checkRangeDouble("TRANS",
                "LongOrig", map_orig_long[n_proj], 1, -180.0, 1, 180.0) != 0)
            ierr = -1;
        if (checkRangeDouble("TRANS",
                "RotCW", map_rot[n_proj], 1, -360.0, 1, 360.0) != 0)
            ierr = -1;
        if (ierr < 0 || istat != 4) {
            nll_puterr("ERROR: reading SIMPLE transformation parameters");
            return (-1);
        }


    } else if (strcmp(map_trans_type[n_proj], "SDC") == 0) {

        map_itype[n_proj] = MAP_TRANS_SDC;
        istat = sscanf(in_line, "%s %lf %lf %lf",
                map_trans_type[n_proj], &map_orig_lat[n_proj], &map_orig_long[n_proj],
                &map_rot[n_proj]);

        angle = -cRPD * map_rot[n_proj];
        map_cosang[n_proj] = cos(angle);
        map_sinang[n_proj] = sin(angle);


        sprintf(MapProjStr[n_proj],
                "TRANSFORM  %s LatOrig %lf  LongOrig %lf  RotCW %lf",
                map_trans_type[n_proj], map_orig_lat[n_proj], map_orig_long[n_proj], map_rot[n_proj]);
        nll_putmsg(3, MapProjStr[n_proj]);

        ierr = 0;
        if (checkRangeDouble("TRANS",
                "LatOrig", map_orig_lat[n_proj], 1, -90.0, 1, 90.0) != 0)
            ierr = -1;
        if (checkRangeDouble("TRANS",
                "LongOrig", map_orig_long[n_proj], 1, -180.0, 1, 180.0) != 0)
            ierr = -1;
        if (checkRangeDouble("TRANS",
                "RotCW", map_rot[n_proj], 1, -360.0, 1, 360.0) != 0)
            ierr = -1;
        if (ierr < 0 || istat != 4) {
            nll_puterr("ERROR: reading SDC transformation parameters");
            return (-1);
        }

        // SDC initialization
        //  conversion factor for latitude
        dlt1 = atan(MAP_TRANS_SDC_DRLT * tan(map_orig_lat[n_proj] * DE2RA));
        dlt2 = atan(MAP_TRANS_SDC_DRLT * tan((map_orig_lat[n_proj] + 1.0) * DE2RA));
        del = dlt2 - dlt1;
        r = ERAD * (1.0 - pow(sin(dlt1), 2) * FLATTENING);
        map_sdc_xltkm[n_proj] = r * del;
        //  conversion factor for longitude
        del = acos(1.0 - (1.0 - cos(DE2RA)) * pow(cos(dlt1), 2));
        bc = r * del;
        map_sdc_xlnkm[n_proj] = bc / cos(dlt1);



    } else if (strcmp(map_trans_type[n_proj], "LAMBERT") == 0) {

        map_itype[n_proj] = MAP_TRANS_LAMBERT;
        istat = sscanf(in_line, "%s %s %lf %lf %lf %lf %lf",
                map_trans_type[n_proj], map_ref_ellipsoid[n_proj],
                &map_orig_lat[n_proj], &map_orig_long[n_proj],
                &map_lambert_1st_std_paral[n_proj], &map_lambert_2nd_std_paral[n_proj],
                &map_rot[n_proj]);

        ierr = 0;
        if (checkRangeDouble("TRANS",
                "LatOrig", map_orig_lat[n_proj], 1, -90.0, 1, 90.0) != 0)
            ierr = -1;
        if (checkRangeDouble("TRANS",
                "LongOrig", map_orig_long[n_proj], 1, -180.0, 1, 180.0) != 0)
            ierr = -1;
        if (checkRangeDouble("TRANS",
                "FirstStdParal", map_lambert_1st_std_paral[n_proj], 1, -90.0, 1, 90.0) != 0)
            ierr = -1;
        if (checkRangeDouble("TRANS",
                "SecondStdParal", map_lambert_2nd_std_paral[n_proj], 1, -90.0, 1, 90.0) != 0)
            ierr = -1;
        if (checkRangeDouble("TRANS",
                "RotCW", map_rot[n_proj], 1, -360.0, 1, 360.0) != 0)
            ierr = -1;

        angle = -cRPD * map_rot[n_proj];
        map_cosang[n_proj] = cos(angle);
        map_sinang[n_proj] = sin(angle);

        /* initialize GMT projection values */
        if (map_setup_proxy(n_proj, map_ref_ellipsoid[n_proj]) < 0) {
            nll_puterr(
                    "ERROR: initializing general transformation parameters, RefEllipsoid may be invalid");
            return (-1);
        }

        /* initialize lambert projection */
        vlamb(n_proj, map_orig_long[n_proj], map_orig_lat[n_proj],
                map_lambert_1st_std_paral[n_proj], map_lambert_2nd_std_paral[n_proj]);

        sprintf(MapProjStr[n_proj],
                "TRANSFORM  %s RefEllipsoid %s  LatOrig %lf  LongOrig %lf  FirstStdParal %lf  SecondStdParal %lf  RotCW %lf",
                map_trans_type[n_proj], map_ref_ellipsoid[n_proj],
                map_orig_lat[n_proj], map_orig_long[n_proj],
                map_lambert_1st_std_paral[n_proj], map_lambert_2nd_std_paral[n_proj],
                map_rot[n_proj]);
        nll_putmsg(3, MapProjStr[n_proj]);

        if (ierr < 0 || istat != 7) {
            nll_puterr("ERROR: reading LAMBERT transformation parameters");
            return (-1);
        }

    } else {

        nll_puterr("ERROR: unrecognized map transformation type");
        return (-1);
    }

    return (0);
}


/** function to convert lat/long to rectangular km coord */
#ifdef _MSC_VER
__declspec(dllexport)
#endif
int latlon2rect(int n_proj, double dlat, double dlong, double* pxrect, double* pyrect) {

    double xtemp, ytemp;

    double xlt1;


    if (map_itype[n_proj] == MAP_TRANS_GLOBAL) {
        *pxrect = dlong;
        *pyrect = dlat;
        return (0);

    } else if (map_itype[n_proj] == MAP_TRANS_NONE) {
        *pxrect = dlong;
        *pyrect = dlat;
        return (0);

    } else if (map_itype[n_proj] == MAP_TRANS_SIMPLE) {
        xtemp = dlong - map_orig_long[n_proj];
        if (xtemp > 180.0)
            xtemp -= 360.0;
        if (xtemp < -180.0)
            xtemp += 360.0;
        xtemp = xtemp * c111 * cos(cRPD * dlat);
        ytemp = (dlat - map_orig_lat[n_proj]) * c111;
        *pxrect = xtemp * map_cosang[n_proj] - ytemp * map_sinang[n_proj];
        *pyrect = ytemp * map_cosang[n_proj] + xtemp * map_sinang[n_proj];
        return (0);

    } else if (map_itype[n_proj] == MAP_TRANS_SDC) {

        /*
        c  now convert lat and lon differences to km
        x=xln-pln
        c     x=pln-xln
        y=plt-xlt
        xlt1=atan(drlt*tan(drad*(plt+xlt)/120.))
        x=x*xlnkm*cos(xlt1)
        y=y*xltkm
         */

        xtemp = dlong - map_orig_long[n_proj];
        if (xtemp > 180.0)
            xtemp -= 360.0;
        if (xtemp < -180.0)
            xtemp += 360.0;
        ytemp = dlat - map_orig_lat[n_proj];

        xlt1 = atan(MAP_TRANS_SDC_DRLT * tan(DE2RA * (dlat + map_orig_lat[n_proj]) / 2.0));
        xtemp = xtemp * map_sdc_xlnkm[n_proj] * cos(xlt1);
        ytemp = ytemp * map_sdc_xltkm[n_proj];

        *pxrect = xtemp * map_cosang[n_proj] - ytemp * map_sinang[n_proj];
        *pyrect = ytemp * map_cosang[n_proj] + xtemp * map_sinang[n_proj];
        return (0);

    } else if (map_itype[n_proj] == MAP_TRANS_LAMBERT) {

        lamb(n_proj, dlong, dlat, &xtemp, &ytemp);
        xtemp /= 1000.0; /* m -> km */
        ytemp /= 1000.0; /* m -> km */
        *pxrect = xtemp * map_cosang[n_proj] - ytemp * map_sinang[n_proj];
        *pyrect = ytemp * map_cosang[n_proj] + xtemp * map_sinang[n_proj];

        return (0);
    }

    return (-1);

}


/** function to convert rectangular km coord to lat/long */
#ifdef _MSC_VER
__declspec(dllexport)
#endif
int rect2latlon(int n_proj, double xrect, double yrect, double* pdlat, double* pdlong) {

    double xtemp, ytemp;

    double xlt1;

    if (map_itype[n_proj] == MAP_TRANS_GLOBAL) {
        *pdlat = yrect;
        *pdlong = xrect;

        return (0);

    } else if (map_itype[n_proj] == MAP_TRANS_NONE) {
        *pdlat = yrect;
        *pdlong = xrect;

        return (0);

    } else if (map_itype[n_proj] == MAP_TRANS_SIMPLE) {
        xtemp = xrect * map_cosang[n_proj] + yrect * map_sinang[n_proj];
        ytemp = yrect * map_cosang[n_proj] - xrect * map_sinang[n_proj];
        *pdlat = map_orig_lat[n_proj] + ytemp / c111;
        *pdlong = map_orig_long[n_proj] + xtemp / (c111 * cos(cRPD * *pdlat));
        // 20121005 AJL - prevent longitude outside of -180 -> 180 deg range
        if (*pdlong < -180.0)
            *pdlong += 360.0;
        else if (*pdlong > 180.0)
            *pdlong -= 360.0;

        return (0);

    } else if (map_itype[n_proj] == MAP_TRANS_SDC) {
        xtemp = xrect * map_cosang[n_proj] + yrect * map_sinang[n_proj];
        ytemp = yrect * map_cosang[n_proj] - xrect * map_sinang[n_proj];

        /*
        c
        fy=fy/xltkm
        plt=xlt+fy
        c
        xlt1=atan(rlt*tan(rad*(plt+xlt)/120.))
        fx=fx/(xlnkm*cos(xlt1))
        pln=xln+fx
         */

        ytemp = ytemp / map_sdc_xltkm[n_proj];
        *pdlat = map_orig_lat[n_proj] + ytemp;
        xlt1 = atan(MAP_TRANS_SDC_DRLT * tan(DE2RA * (*pdlat + map_orig_lat[n_proj]) / 2.0));
        xtemp = xtemp / (map_sdc_xlnkm[n_proj] * cos(xlt1));
        *pdlong = map_orig_long[n_proj] + xtemp;
        // 20121005 AJL - prevent longitude outside of -180 -> 180 deg range
        if (*pdlong < -180.0)
            *pdlong += 360.0;
        else if (*pdlong > 180.0)
            *pdlong -= 360.0;

        return (0);

    } else if (map_itype[n_proj] == MAP_TRANS_LAMBERT) {

        xtemp = xrect * map_cosang[n_proj] + yrect * map_sinang[n_proj];
        ytemp = yrect * map_cosang[n_proj] - xrect * map_sinang[n_proj];
        ilamb(n_proj, pdlong, pdlat, xtemp * 1000.0, ytemp * 1000.0);
        // 20121005 AJL - prevent longitude outside of -180 -> 180 deg range
        if (*pdlong < -180.0)
            *pdlong += 360.0;
        else if (*pdlong > 180.0)
            *pdlong -= 360.0;

        return (0);
    }

    return (-1);
}


/** function to convert rectangular km angle to lat/long angle */
double rect2latlonAngle(int n_proj, double rectAngle) {
    double angle;

    if (map_itype[n_proj] == MAP_TRANS_SIMPLE ||
            map_itype[n_proj] == MAP_TRANS_SDC ||
            map_itype[n_proj] == MAP_TRANS_LAMBERT) {
        angle = rectAngle - map_rot[n_proj];
        if (angle < 0.0)
            angle += 360.0;
        else if (angle > 360.0)
            angle -= 360.0;
        return (angle);
    } else
        return (rectAngle);
}


/** function to convert lat/long km angle to rectangular angle */
double latlon2rectAngle(int n_proj, double latlonAngle) {
    double angle;

    if (map_itype[n_proj] == MAP_TRANS_SIMPLE ||
            map_itype[n_proj] == MAP_TRANS_SDC ||
            map_itype[n_proj] == MAP_TRANS_LAMBERT) {
        angle = latlonAngle + map_rot[n_proj];
        if (angle < 0.0)
            angle += 360.0;
        else if (angle > 360.0)
            angle -= 360.0;
        return (angle);
    } else
        return (latlonAngle);
}
