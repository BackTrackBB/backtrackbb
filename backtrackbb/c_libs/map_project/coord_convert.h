#ifndef COORD_CONVERT_H
#define COORD_CONVERT_H

#include <math.h>
#include <string.h>
#include "geo.h"
#include "map_project.h"
#include "util.h"

#ifndef MAXLINE
#define MAXLINE 101
#endif

// mode
#define MODE_RECT 0 // rectangular cartesian x(km),y(km),z:depth(km)
#define MODE_GLOBAL 1 // spherical x:longitdue(deg),y:latittude(deg),z:depth(km)
int GeometryMode;

/* geographic transformations (lat/long <=> x/y) */
#define NUM_PROJ_MAX 3
#define MAP_TRANS_UNDEF -1
#define MAP_TRANS_NONE 0
#define MAP_TRANS_GLOBAL 1
#define MAP_TRANS_SIMPLE 2
#define MAP_TRANS_LAMBERT 3
#define MAP_TRANS_SDC 4
char map_trans_type[NUM_PROJ_MAX][MAXLINE]; /* name of projection */
int map_itype[NUM_PROJ_MAX]; /* int id of projection */
char MapProjStr[NUM_PROJ_MAX][2 * MAXLINE]; /* string description of proj params */
char map_ref_ellipsoid[NUM_PROJ_MAX][MAXLINE]; /* name of reference ellipsoid */
/* general map parameters */
double map_orig_lat[NUM_PROJ_MAX], map_orig_long[NUM_PROJ_MAX], map_rot[NUM_PROJ_MAX];
double map_cosang[NUM_PROJ_MAX], map_sinang[NUM_PROJ_MAX]; /* rotation */
/* LAMBERT projection parameters */
double map_lambert_1st_std_paral[NUM_PROJ_MAX], map_lambert_2nd_std_paral[NUM_PROJ_MAX];
/* SDC Short Distance Coversion projection parameters */
double map_sdc_xltkm[NUM_PROJ_MAX], map_sdc_xlnkm[NUM_PROJ_MAX];
#define MAP_TRANS_SDC_DRLT 0.99330647

/* constants */
double cPI;
double cRPD;
double c111;

/** function to read map transformation parameters from input line ***/
#ifdef _MSC_VER
__declspec(dllexport)
#endif
int get_transform(int n_proj, char* in_line);
/** function to convert lat/long to rectangular km coord */
#ifdef _MSC_VER
__declspec(dllexport)
#endif
int latlon2rect(int n_proj, double dlat, double dlong, double* pxrect, double* pyrect);
/** function to convert rectangular km coord to lat/long */
#ifdef _MSC_VER
__declspec(dllexport)
#endif
int rect2latlon(int n_proj, double xrect, double yrect, double* pdlat, double* pdlong);
/** function to convert rectangular km angle to lat/long angle */
double rect2latlonAngle(int n_proj, double rectAngle);
/** function to convert lat/long km angle to rectangular angle */
double latlon2rectAngle(int n_proj, double latlonAngle);

#endif
