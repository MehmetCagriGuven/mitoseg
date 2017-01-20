/*
 * curve.h
 *
 *  Created on: Jan 13, 2017
 *      Author: fst
 */

#ifndef CURVE_H_
#define CURVE_H_

#include "opencv/cv.h"
#include "opencv/highgui.h"

#include "preprocess.h"

#include <stdio.h>

typedef struct {
	float energy[LE_BINS];
	int major;
} localEnergy;

typedef struct {
	// Geometric model params
	int x1, y1, x2, y2, h;

	// Algebraic model params
	double R, a, b, sint, cost;

	// Curve fitness params
	float score, avgscore;
	int len;

} curve;

void getLocalEnergyMap(ridge *r, const int width, const int height,
		const int ssize, const int bsize, localEnergy **e, int *w, int *h);
void visualizeLEMajority(localEnergy *e, IplImage **visual, const int w,
		const int h, const int ssize, const int showDir);
void calcCurveParams(localEnergy *e, const int w, const int h, curve *c,
		const double S);
void calcAuxCurveParams(curve *c);
void drawCurve(IplImage *im, curve *c, const int ssize, const int psize,
		const double S, CvScalar cl = cvScalar(255, 0, 0));
void drawCurve2(IplImage *im, curve *c, const int ssize, const int psize,
		const double S, CvScalar cl = cvScalar(255, 0, 0));
int fitCurve3(localEnergy *e, const int w, const int h, curve *c,
		const int maxIter, char *mask = NULL);
int fitCurve2(localEnergy *e, const int w, const int h, curve *c,
		const int maxIter, char *mask = NULL);
int fitCurve(localEnergy *e, const int w, const int h, curve *c,
		const int maxIter, char *mask = NULL);
int fitCurves(localEnergy *e, const int w, const int h, curve **clist,
		const int FCS_MAXCURVE, const int FCS_XYSTEP, const int FCS_XYRANGE,
		const float FCS_INIT_THRESH, const int FCS_MINLEN);
void visualizeCurves(IplImage *bckgnd, curve *clist, const int nc,
		const int ssize, IplImage *visual, const float FLC_THRESH = -1.0f,
		const float FLC_AVGTHRESH = -1.0f);
void initCurveFitting();
void drawBinaryCurves(curve *clist, const int nc, const int psize,
		IplImage **out, const int w, const int h);
void visualizeHiLoCurves(curve *clist_lo, const int nc_lo, const int w_lo,
		const int h_lo, curve *clist_hi, const int nc_hi,
		const double maxoverlap, IplImage *bckgnd, IplImage *visual);
void getHixCurves(curve *clist_lo, const int nc_lo, const int w_lo,
		const int h_lo, curve *clist_hi, const int nc_hi, curve **clist_hilo,
		int *nc_hilo);
void filterCurves(const curve *cin, const int ncin, curve **cout, int *ncout,
		const float FLC_THRESH = -1.0f, const float FLC_AVGTHRESH = -1.0f);
void getCurvePoints(const curve *c, const double S, double *cx, double *cy);
void saveLocalEnergyStructure(int s, char const *tag, localEnergy *e, int w,
		int h);
void loadLocalEnergyStructure(int s, char const *tag, localEnergy **e, int *w,
		int *h);
void saveCurveStructure(int s, char const *tag, curve *c, int numRecords);
int loadCurveStructure(int s, char const *tag, curve **c);
int isInside(const curve *c, double x, double y, double bias = 0.0);

#endif /* CURVE_H_ */
