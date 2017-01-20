/*
 * snake25d.h
 *
 *  Created on: Jan 15, 2017
 *      Author: fst
 */

#ifndef SNAKE25D_H_
#define SNAKE25D_H_

#include "opencv/cv.h"
#include "opencv/highgui.h"

#include "curve.h"

#define SNAKE_GROWING 0
#define SNAKE_CONVERGED 1
#define SNAKE_MAXITER 2
#define SNAKE_MAXAREA 3
#define SNAKE_MINAREA 4
//
#define SNAKE_INVALID 0
#define SNAKE_VALID 1
//

struct snake25d {
	double node[SN25D_T][SN_N][2];
	double d1_xy[SN25D_T][SN_N][2];
	double d2_xy[SN25D_T][SN_N][2];
	double d4_xy[SN25D_T][SN_N][2];
	double d1_z[SN25D_T][SN_N][2];
	double d2_z[SN25D_T][SN_N][2];
	double d4_z[SN25D_T][SN_N][2];
	double grad_Eint[SN25D_T][SN_N][2];
	double grad_Eext[SN25D_T][SN_N][2];
	double area[SN25D_T];
	char status;
	char isValid[SN25D_T];
	double validity;
	int start_z, end_z, t;
	double aspectRatio_z;
};

struct snake25dList {
	snake25d snake;
	snake25dList *next;
};

struct curveList {
	curve *clist;
	int n;
};

struct curveStack {
	curveList **cstack;
	int start_z, end_z, t;
};

template<typename T>
struct curvePoints {
	T *cx;
	T *cy;
	int len;
};

template<typename T>
struct curvePointsStack {
	curvePoints<T> **cps;
	int start_z, end_z, t;
};

struct curveEnergyMap {
	double *ce;
	int w, h;
};

struct curveEnergyStack {
	curveEnergyMap **ces;
	int start_z, end_z, t;
};

struct curveEnergyGradientMap {
	double (*grad_ce)[2];
	int w, h;
};

struct curveEnergyGradientStack {
	curveEnergyGradientMap **grad_ces;
	int start_z, end_z, t;
};

struct dataPacket {
	curveStack *cs_lo;
	curveStack *cs_hi;
	curveEnergyStack *ces_lo;
	curveEnergyStack *ces_hi;
	curveEnergyGradientStack *grad_ces_lo;
	curvePointsStack<double> *cps_lo;
	curvePointsStack<int> *cps_hi;
};
//
curveStack* loadCurveStack(int start_z, int end_z, char const *tag);
void destroyCurveStack(curveStack *cstack);
curvePoints<double>* createCurvePoints(curveList *c);
curvePoints<int>* createCurvePoints(curveList *c, int w, int h);
curvePointsStack<double>* createCurvePointsStack(curveStack *cstack);
curvePointsStack<int>* createCurvePointsStack(curveStack *cstack, int w, int h);
void destroyCurvePointsStack(curvePointsStack<int> *cps);
curveEnergyStack* createCurveEnergyStack(curveStack *cs, const int w,
		const int h);
void destroyCurveEnergyStack(curveEnergyStack *ces);
curveEnergyGradientStack* createCurveEnergyGradientStack(curveEnergyStack *ces);
void destroyCurveEnergyGradientStack(curveEnergyGradientStack *grad_ces);
//
void getMapSize(CvSize ridgeImageSize, const int ssize, int *w, int *h);
void initSnake25d(snake25d *s, const double x, const double y,
		const int start_z, const int end_z, const double init_r, const int w,
		const int h, const double ar_z);
void addSnake(snake25dList **slist, snake25d *snake);
void removeSnake(snake25dList **slist);
void destroySnakeList(snake25dList **slist);
int convertSnakeList2Array(snake25dList *slist, snake25d **sarray);
void startSnake25d(curveEnergyGradientStack *grad_ces, snake25d *outputSnake,
		const double w_tens_xy, const double w_curv_xy, const double w_tens_z,
		const double w_curv_z, const double w_curve, const double w_inf,
		const double aspectRatio_z, snake25d *initialSnake, const double initx,
		const double inity, const int start_z, int end_z);
int inflationSnake25d(dataPacket *dp, const double initx, const double inity,
		const int start_z, int end_z, const double w_tens_xy,
		const double w_curv_xy, const double w_tens_z, const double w_curv_z,
		const double w_curve, const double w_inf_min, const double w_inf_max,
		const double w_inf_step, const double aspectRatio_z,
		snake25dList **outputSnakeList);
void reportValidation(dataPacket *dp, snake25d *s);
void validateSnake25d(dataPacket *dp, snake25d *s);
void validateSnakeList(dataPacket *dp, snake25dList *slist);
int findInitialPointsWithinZRange(dataPacket *dp, CvPoint **pts, int start_z,
		int end_z);
int retrieveAllSnakes25d(dataPacket *dp, CvPoint *initPts, int n, int start_z,
		int end_z, snake25dList **outputSnakeList);
void saveSnakeArray(snake25d *sarray, int n, int tag_start_z, int tag_end_z);
int loadSnakeArray(int tag_start_z, int tag_end_z, snake25d **sarray);
void savePoints(CvPoint *pts, int n, int tag_start_z, int tag_end_z);
int loadPoints(int tag_start_z, int tag_end_z, CvPoint **pts);
int filterSnakeArrayByValidity(snake25d *in, int n, double th, snake25d **out);
void filterSnake25dByValidity(snake25dList *in, snake25dList **out, double th);
int disintegrateSnake25dTo2d(snake25d *in, snake25d *outArray);
int disintegrateArrayOfSnake25dTo2d(snake25d *sarray, int n,
		snake25d **outArray);
int filterArrayOfSnake25dBySimilarity(snake25d *sarray, int n,
		snake25d **outArray, double th);
void getEnergyGradient(const double *ce, const int cw, const int ch,
		double (*grad_ce)[2]);
void getCurveEnergyImage(curve *clist, const int n, const int w, const int h,
		double **ce);
#endif /* SNAKE25D_H_ */
