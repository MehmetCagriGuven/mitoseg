/*
 * snake25d.cpp
 *
 *  Created on: Jan 15, 2017
 *      Author: fst
 */

#include "snake25d.h"
#include "detection.h"
#include <pthread.h>
#include <semaphore.h>
#include <math.h>

sem_t sem2;
pthread_mutex_t addlock = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t tnslock = PTHREAD_MUTEX_INITIALIZER;

curveList* loadCurveList(int s, char const *tag) {
	curveList *cl = (curveList *) malloc(sizeof(curveList));
	cl->clist = NULL;
	cl->n = loadCurveStructure(s, tag, &cl->clist);
	return cl;
}

void destroyCurveList(curveList *cl) {
	if (cl) {
		if (cl->clist)
			free(cl->clist);
		free(cl);
	}
}

curveStack* loadCurveStack(int start_z, int end_z, char const *tag) {
	int s, i;
	curveStack *cstack = (curveStack *) malloc(sizeof(curveStack));
	cstack->start_z = start_z;
	cstack->end_z = end_z;
	cstack->t = end_z - start_z + 1;
	cstack->cstack = (curveList **) malloc(sizeof(curveList *) * cstack->t);

	for (i = 0, s = start_z; s <= end_z; i++, s++) {
		cstack->cstack[i] = loadCurveList(s, tag);
	}
	return cstack;
}

void destroyCurveStack(curveStack *cstack) {
	int s, i;
	if (cstack) {
		for (i = 0, s = cstack->start_z; s <= cstack->end_z; i++, s++) {
			destroyCurveList(cstack->cstack[i]);
		}
		free(cstack);
	}
}

curvePoints<double>* createCurvePoints(curveList *c) {
	int nc = 0;
	int i, t;

	for (i = 0; i < c->n; i++)
		nc += c->clist[i].len;

	curvePoints<double> *cp = (curvePoints<double> *) malloc(
			sizeof(curvePoints<double> ));
	cp->cx = (double *) malloc(sizeof(double) * nc);
	cp->cy = (double *) malloc(sizeof(double) * nc);
	cp->len = nc;

	for (i = 0, t = 0; i < c->n; t += c->clist[i++].len) {
		getCurvePoints(&c->clist[i], 1.0, &cp->cx[t], &cp->cy[t]);
	}
	return cp;
}

template<typename T>
void destroyCurvePoints(curvePoints<T> *cp) {
	if (cp) {
		if (cp->cx)
			free(cp->cx);
		if (cp->cy)
			free(cp->cy);
		free(cp);
	}
}

curvePoints<int>* createCurvePoints(curveList *c, int w, int h) {
	curvePoints<double> *temp = createCurvePoints(c);

	curvePoints<int> *cp = (curvePoints<int> *) malloc(
			sizeof(curvePoints<int> ));
	cp->cx = (int *) malloc(sizeof(int) * temp->len);
	cp->cy = (int *) malloc(sizeof(int) * temp->len);
	cp->len = temp->len;

	int j;
	int w_1 = w - 1;
	int h_1 = h - 1;
	for (j = 0; j < cp->len; j++) {
		cp->cx[j] = cvRound(temp->cx[j]);
		cp->cy[j] = cvRound(temp->cy[j]);
		if (cp->cx[j] < 0)
			cp->cx[j] = 0;
		if (cp->cy[j] < 0)
			cp->cy[j] = 0;
		if (cp->cx[j] > w_1)
			cp->cx[j] = w_1;
		if (cp->cy[j] > h_1)
			cp->cy[j] = h_1;
	}
	destroyCurvePoints(temp);

	int msize = w * h;
	char *map = (char *) malloc(sizeof(char) * msize);
	memset(map, 0, sizeof(char) * msize);
	for (j = 0; j < cp->len; j++)
		map[cp->cy[j] * w + cp->cx[j]] = 1;
	cp->len = 0;
	for (j = 0; j < msize; j++) {
		if (map[j]) {
			cp->cx[cp->len] = j % w;
			cp->cy[cp->len] = j / w;
			cp->len++;
		}
	}
	free(map);

	return cp;
}

curvePointsStack<double>* createCurvePointsStack(curveStack *cstack) {
	curvePointsStack<double> *cps = (curvePointsStack<double> *) malloc(
			sizeof(curvePointsStack<double> ));
	cps->cps = (curvePoints<double> **) malloc(
			sizeof(curvePoints<double> *) * cstack->t);
	cps->start_z = cstack->start_z;
	cps->end_z = cstack->end_z;
	cps->t = cstack->t;
	for (int i = 0; i < cps->t; i++) {
		cps->cps[i] = createCurvePoints(cstack->cstack[i]);
	}
	return cps;
}

curvePointsStack<int>* createCurvePointsStack(curveStack *cstack, int w,
		int h) {
	curvePointsStack<int> *cps = (curvePointsStack<int> *) malloc(
			sizeof(curvePointsStack<int> ));
	cps->cps = (curvePoints<int> **) malloc(
			sizeof(curvePoints<int> *) * cstack->t);
	cps->start_z = cstack->start_z;
	cps->end_z = cstack->end_z;
	cps->t = cstack->t;
	for (int i = 0; i < cps->t; i++) {
		cps->cps[i] = createCurvePoints(cstack->cstack[i], w, h);
	}
	return cps;
}

template<typename T>
void _destroyCurvePointsStack(curvePointsStack<T> *cps) {
	if (cps) {
		for (int i = 0; i < cps->t; i++) {
			destroyCurvePoints(cps->cps[i]);
		}
		free(cps);
	}
}

void destroyCurvePointsStack(curvePointsStack<int> *cps) {
	_destroyCurvePointsStack(cps);
}

curveEnergyMap* createCurveEnergyMap(curveList *cl, const int w, const int h) {
	curveEnergyMap *ce = (curveEnergyMap *) malloc(sizeof(curveEnergyMap));
	ce->w = w;
	ce->h = h;
	ce->ce = NULL;
	getCurveEnergyImage(cl->clist, cl->n, w, h, &ce->ce);
	return ce;
}

void destroyCurveEnergyMap(curveEnergyMap *ce) {
	if (ce) {
		if (ce->ce)
			free(ce->ce);
		free(ce);
	}
}

curveEnergyStack* createCurveEnergyStack(curveStack *cs, const int w,
		const int h) {
	curveEnergyStack *ces = (curveEnergyStack *) malloc(
			sizeof(curveEnergyStack));
	ces->ces = (curveEnergyMap **) malloc(sizeof(curveEnergyMap *) * cs->t);
	ces->start_z = cs->start_z;
	ces->end_z = cs->end_z;
	ces->t = cs->t;

	for (int i = 0; i < ces->t; i++) {
		ces->ces[i] = createCurveEnergyMap(cs->cstack[i], w, h);
	}
	return ces;
}

void destroyCurveEnergyStack(curveEnergyStack *ces) {
	if (ces) {
		for (int i = 0; i < ces->t; i++) {
			if (ces->ces[i]) {
				destroyCurveEnergyMap(ces->ces[i]);
			}
		}
		free(ces);
	}
}

curveEnergyGradientMap* createCurveEnergyGradientMap(curveEnergyMap *ce) {
	curveEnergyGradientMap *grad_ce = (curveEnergyGradientMap *) malloc(
			sizeof(curveEnergyGradientMap));
	int msize = ce->w * ce->h;
	grad_ce->grad_ce = (double (*)[2]) malloc(sizeof(double[2]) * msize);
	grad_ce->w = ce->w;
	grad_ce->h = ce->h;
	getEnergyGradient(ce->ce, ce->w, ce->h, grad_ce->grad_ce);
	return grad_ce;
}

void destroyCurveEnergyGradientMap(curveEnergyGradientMap *grad_ce) {
	if (grad_ce) {
		if (grad_ce->grad_ce)
			free(grad_ce->grad_ce);
		free(grad_ce);
	}
}

curveEnergyGradientStack* createCurveEnergyGradientStack(
		curveEnergyStack *ces) {
	curveEnergyGradientStack *grad_ces = (curveEnergyGradientStack *) malloc(
			sizeof(curveEnergyGradientStack));
	grad_ces->grad_ces = (curveEnergyGradientMap **) malloc(
			sizeof(curveEnergyGradientMap *) * ces->t);
	grad_ces->start_z = ces->start_z;
	grad_ces->end_z = ces->end_z;
	grad_ces->t = ces->t;

	for (int i = 0; i < ces->t; i++) {
		grad_ces->grad_ces[i] = createCurveEnergyGradientMap(ces->ces[i]);
	}
	return grad_ces;
}

void destroyCurveEnergyGradientStack(curveEnergyGradientStack *grad_ces) {
	if (grad_ces) {
		for (int i = 0; i < grad_ces->t; i++) {
			if (grad_ces->grad_ces[i]) {
				destroyCurveEnergyGradientMap(grad_ces->grad_ces[i]);
			}
		}
		free(grad_ces);
	}
}

void getMapSize(CvSize ridgeImageSize, const int ssize, int *w, int *h) {
	*w = (ridgeImageSize.width / ssize)
			+ ((ridgeImageSize.width % ssize) ? 1 : 0);
	*h = (ridgeImageSize.height / ssize)
			+ ((ridgeImageSize.height % ssize) ? 1 : 0);
}

void compSnakeSliceArea(snake25d *s, const int z) {
	int k, j;
	double area = 0;
	for (k = 0, j = 1; k < SN_N; k++, j = (j + 1) % SN_N) {
		area += s->node[z][k][0] * s->node[z][j][1]
				- s->node[z][k][1] * s->node[z][j][0];
	}
	s->area[z] = 0.5 * area;
	if (s->area[z] < 0)
		s->area[z] = -s->area[z];
}

void compSnakeArea(snake25d *s) {
	for (int z = 0; z < s->t; z++) {
		compSnakeSliceArea(s, z);
	}
}

void compApproxIntersectionArea(double p1[][2], int n1, double p2[][2], int n2,
		double *a1, double *a2, double *is, double precision = 1.0) {
	int i, x1, y1, x2, y2, w, h;
	CvPoint *k1 = NULL;
	CvPoint *k2 = NULL;
	if (n1)
		k1 = (CvPoint *) malloc(sizeof(CvPoint) * n1);
	if (n2)
		k2 = (CvPoint *) malloc(sizeof(CvPoint) * n2);
	for (i = 0; i < n1; i++) {
		k1[i].x = (int) (p1[i][0] * precision);
		k1[i].y = (int) (p1[i][1] * precision);
	}
	for (i = 0; i < n2; i++) {
		k2[i].x = (int) (p2[i][0] * precision);
		k2[i].y = (int) (p2[i][1] * precision);
	}
	if (n1) {
		x1 = x2 = k1[0].x;
		y1 = y2 = k1[0].y;
	} else {
		x1 = x2 = k2[0].x;
		y1 = y2 = k2[0].y;
	}
	for (i = 0; i < n1; i++) {
		x1 = (x1 > k1[i].x) ? k1[i].x : x1;
		y1 = (y1 > k1[i].y) ? k1[i].y : y1;
		x2 = (x2 < k1[i].x) ? k1[i].x : x2;
		y2 = (y2 < k1[i].y) ? k1[i].y : y2;
	}
	for (i = 0; i < n2; i++) {
		x1 = (x1 > k2[i].x) ? k2[i].x : x1;
		y1 = (y1 > k2[i].y) ? k2[i].y : y1;
		x2 = (x2 < k2[i].x) ? k2[i].x : x2;
		y2 = (y2 < k2[i].y) ? k2[i].y : y2;
	}
	w = x2 - x1 + 3;
	h = y2 - y1 + 3;
	for (i = 0; i < n1; i++) {
		k1[i].x -= (x1 - 1);
		k1[i].y -= (y1 - 1);
	}
	for (i = 0; i < n2; i++) {
		k2[i].x -= (x1 - 1);
		k2[i].y -= (y1 - 1);
	}
	IplImage *buf1 = cvCreateImage(cvSize(w, h), IPL_DEPTH_8U, 1);
	IplImage *buf2 = cvCreateImage(cvSize(w, h), IPL_DEPTH_8U, 1);
	cvSetZero(buf1);
	cvSetZero(buf2);
	if (n1)
		cvFillPoly(buf1, &k1, &n1, 1, cvScalar(255, 255, 255));
	if (n2)
		cvFillPoly(buf2, &k2, &n2, 1, cvScalar(255, 255, 255));
	int c1 = cvCountNonZero(buf1);
	int c2 = cvCountNonZero(buf2);
	cvAnd(buf1, buf2, buf1);
	int inter = cvCountNonZero(buf1);
	*a1 = c1 / (precision * precision);
	*a2 = c2 / (precision * precision);
	*is = inter / (precision * precision);
	if (k1)
		free(k1);
	if (k2)
		free(k2);
	cvReleaseImage(&buf1);
	cvReleaseImage(&buf2);
}

double getDiceSimilarity2d(double p1[][2], int n1, double p2[][2], int n2,
		double precision = 1.0) {
	double a1, a2, is, ds;
	compApproxIntersectionArea(p1, n1, p2, n2, &a1, &a2, &is, precision);
	ds = 2.0 * is / (a1 + a2);
	return ds;
}

double getDiceSimilaritySnake25d(snake25d *s1, snake25d *s2, double precision =
		1.0) {
	int z1 = (s1->start_z < s2->start_z) ? s1->start_z : s2->start_z;
	int z2 = (s1->end_z > s2->end_z) ? s1->end_z : s2->end_z;
	int z;
	double a1, a2, is;
	double v1 = 0, v2 = 0, iv = 0, ds = 1.0;

	for (z = z1; z <= z2; z++) {
		if (z >= s1->start_z && z >= s2->start_z && z <= s1->end_z
				&& z <= s2->end_z) {
			compApproxIntersectionArea(s1->node[z - s1->start_z], SN_N,
					s2->node[z - s2->start_z], SN_N, &a1, &a2, &is, precision);
		} else if (z >= s1->start_z && z <= s1->end_z) {
			compApproxIntersectionArea(s1->node[z - s1->start_z], SN_N, NULL, 0,
					&a1, &a2, &is, precision);
		} else {
			compApproxIntersectionArea(NULL, 0, s2->node[z - s2->start_z], SN_N,
					&a1, &a2, &is, precision);
		}
		v1 += a1;
		v2 += a2;
		iv += is;
	}
	ds = 2.0 * iv / (v1 + v2);
	return ds;
}

void initSnake25d(snake25d *s, const double x, const double y,
		const int start_z, const int end_z, const double init_r, const int w,
		const int h, const double ar_z) {
	int i, z;
	double rad, sx, sy;

	s->start_z = start_z;
	s->end_z = end_z;
	s->t = end_z - start_z + 1;
	s->aspectRatio_z = ar_z;
	s->status = SNAKE_GROWING;
	s->validity = 0;

	for (i = 0; i < SN_N; i++) {
		rad = 2 * PI * i / SN_N;
		sx = x + init_r * cos(rad);
		sy = y + init_r * sin(rad);
		// Viewpoint constraints & correction
		if (sx < 0)
			sx = 0;
		if (sy < 0)
			sy = 0;
		if (sx >= w)
			sx = w - 1;
		if (sy >= h)
			sy = h - 1;

		for (z = 0; z < s->t; z++) {
			s->node[z][i][0] = sx;
			s->node[z][i][1] = sy;
			//
			s->grad_Eint[z][i][0] = 0;
			s->grad_Eint[z][i][1] = 0;
			s->grad_Eext[z][i][0] = 0;
			s->grad_Eext[z][i][1] = 0;
			s->d1_xy[z][i][0] = 0;
			s->d1_xy[z][i][1] = 0;
			s->d2_xy[z][i][0] = 0;
			s->d2_xy[z][i][1] = 0;
			s->d4_xy[z][i][0] = 0;
			s->d4_xy[z][i][1] = 0;
			s->d1_z[z][i][0] = 0;
			s->d1_z[z][i][1] = 0;
			s->d2_z[z][i][0] = 0;
			s->d2_z[z][i][1] = 0;
			s->d4_z[z][i][0] = 0;
			s->d4_z[z][i][1] = 0;
			//
			s->isValid[z] = SNAKE_INVALID;
		}
	}

	// Init area
	compSnakeArea(s);
}

void compSnakeDerivatives(snake25d *s) {
	int i, t;
	int k1, k2, k3, k4, k5;
	int t1, t2, t3, t4, t5;

	for (t = 0; t < s->t; t++) {
		// z indices
		t1 = t - 2;
		t2 = t - 1;
		t3 = t;
		t4 = t + 1;
		t5 = t + 2;
		if (t1 < 0)
			t1 = 0;
		if (t2 < 0)
			t2 = 0;
		if (t4 >= s->t)
			t4 = s->t - 1;
		if (t5 >= s->t)
			t5 = s->t - 1;

		// init xy indices
		k1 = SN_N - 2;
		k2 = SN_N - 1;
		k3 = 0;
		k4 = 1;
		k5 = 2;
		for (i = 0; i < SN_N; i++) {
			// xy derivatives
			s->d1_xy[t][i][0] = 0.5 * (s->node[t][k2][0] - s->node[t][k4][0]);
			s->d1_xy[t][i][1] = 0.5 * (s->node[t][k2][1] - s->node[t][k4][1]);
			s->d2_xy[t][i][0] = 0.25
					* (-s->node[t][k2][0] + 2 * s->node[t][k3][0]
							- s->node[t][k4][0]);
			s->d2_xy[t][i][1] = 0.25
					* (-s->node[t][k2][1] + 2 * s->node[t][k3][1]
							- s->node[t][k4][1]);
			s->d4_xy[t][i][0] = 0.0625
					* (s->node[t][k1][0] - 4 * s->node[t][k2][0]
							+ 6 * s->node[t][k3][0] - 4 * s->node[t][k4][0]
							+ s->node[t][k5][0]);
			s->d4_xy[t][i][1] = 0.0625
					* (s->node[t][k1][1] - 4 * s->node[t][k2][1]
							+ 6 * s->node[t][k3][1] - 4 * s->node[t][k4][1]
							+ s->node[t][k5][1]);

			// z derivatives
			s->d1_z[t][i][0] = s->aspectRatio_z * 0.5
					* (s->node[t2][i][0] - s->node[t4][i][0]);
			s->d1_z[t][i][1] = s->aspectRatio_z * 0.5
					* (s->node[t2][i][1] - s->node[t4][i][1]);
			s->d2_z[t][i][0] = s->aspectRatio_z * 0.25
					* (-s->node[t2][i][0] + 2 * s->node[t3][i][0]
							- s->node[t4][i][0]);
			s->d2_z[t][i][1] = s->aspectRatio_z * 0.25
					* (-s->node[t2][i][1] + 2 * s->node[t3][i][1]
							- s->node[t4][i][1]);
			s->d4_z[t][i][0] = s->aspectRatio_z * 0.0625
					* (s->node[t1][i][0] - 4 * s->node[t2][i][0]
							+ 6 * s->node[t3][i][0] - 4 * s->node[t4][i][0]
							+ s->node[t5][i][0]);
			s->d4_z[t][i][1] = s->aspectRatio_z * 0.0625
					* (s->node[t1][i][1] - 4 * s->node[t2][i][1]
							+ 6 * s->node[t3][i][1] - 4 * s->node[t4][i][1]
							+ s->node[t5][i][1]);

			// update xy indices
			k1 = (k1 + 1) % SN_N;
			k2 = (k2 + 1) % SN_N;
			k3 = (k3 + 1) % SN_N;
			k4 = (k4 + 1) % SN_N;
			k5 = (k5 + 1) % SN_N;
		}
	}
}

void compIntEnergyGrad(snake25d *s, const double w_tens_xy,
		const double w_curv_xy, const double w_tens_z, const double w_curv_z) {
	int i, t;
	for (t = 0; t < s->t; t++) {
		for (i = 0; i < SN_N; i++) {
			// Gradient vector of internal energy
			s->grad_Eint[t][i][0] = w_tens_xy * s->d2_xy[t][i][0]
					+ w_curv_xy * s->d4_xy[t][i][0]
					+ w_tens_z * s->d2_z[t][i][0] + w_curv_z * s->d4_z[t][i][0];

			s->grad_Eint[t][i][1] = w_tens_xy * s->d2_xy[t][i][1]
					+ w_curv_xy * s->d4_xy[t][i][1]
					+ w_tens_z * s->d2_z[t][i][1] + w_curv_z * s->d4_z[t][i][1];
		}
	}
}

void compExtEnergyGrad(snake25d *s, curveEnergyGradientStack *grad_ces,
		const double w_curve, const double w_inf) {
	int i, j, x, y, t, k;
	double m, nx, ny;

	for (t = 0, k = s->start_z - grad_ces->start_z; t < s->t; t++, k++) {
		for (i = 0; i < SN_N; i++) {
			x = cvRound(s->node[t][i][0]);
			y = cvRound(s->node[t][i][1]);
			j = y * grad_ces->grad_ces[k]->w + x;
			m = sqrt(
					s->d1_xy[t][i][0] * s->d1_xy[t][i][0]
							+ s->d1_xy[t][i][1] * s->d1_xy[t][i][1]);
			nx = -s->d1_xy[t][i][1] / m;
			ny = s->d1_xy[t][i][0] / m;
			// Compute the gradient of external energy
			s->grad_Eext[t][i][0] = -w_curve
					* grad_ces->grad_ces[k]->grad_ce[j][0] - w_inf * nx;
			s->grad_Eext[t][i][1] = -w_curve
					* grad_ces->grad_ces[k]->grad_ce[j][1] - w_inf * ny;
		}
	}
}

void equidistantCorrectionUV(snake25d *s, double uv[][SN_N][2],
		double uv_corrected[][SN_N][2]) {
	int t, j, jp, jn;
	double dx, dy, dm, cx, cy, c;

	for (t = 0; t < s->t; t++) {
		for (j = 0, jp = SN_N - 1, jn = 1; j < SN_N;
				j++, jp = (jp + 1) % SN_N, jn = (jn + 1) % SN_N) {
			dx = s->node[t][jn][0] - s->node[t][jp][0];
			dy = s->node[t][jn][1] - s->node[t][jp][1];
			dm = dx * dx + dy * dy;
			if (dm < EPS) {
				c = 0;
			} else {
				cx = 0.5 * (s->node[t][jn][0] + s->node[t][jp][0]);
				cy = 0.5 * (s->node[t][jn][1] + s->node[t][jp][1]);
				c = (dx * (cx - (s->node[t][j][0] + uv[t][j][0]))
						+ dy * (cy - (s->node[t][j][1] + uv[t][j][1]))) / dm;
			}
			uv_corrected[t][j][0] = uv[t][j][0] + c * dx;
			uv_corrected[t][j][1] = uv[t][j][1] + c * dy;
		}
	}
}

void compUpdateVectors(snake25d *s, double uv[][SN_N][2]) {
	int t, j;
	double ux, uy, m, mx, gamma;
	for (t = 0; t < s->t; t++) {
		mx = 0;
		for (j = 0; j < SN_N; j++) {
			uv[t][j][0] = ux = s->grad_Eint[t][j][0] + s->grad_Eext[t][j][0];
			uv[t][j][1] = uy = s->grad_Eint[t][j][1] + s->grad_Eext[t][j][1];
			m = ux * ux + uy * uy;
			if (mx < m)
				mx = m;
		}
		mx = sqrt(mx);
		// Compute time discretization parameter gamma
		if (mx > EPS)
			gamma = 1.0 / mx;
		else
			gamma = 0;
		for (j = 0; j < SN_N; j++) {
			uv[t][j][0] *= -SN25D_K * gamma;
			uv[t][j][1] *= -SN25D_K * gamma;
		}
	}
}

void copySnake(snake25d *source, snake25d *dest) {
	memcpy((void *) dest, (void *) source, sizeof(snake25d));
}

void addSnake(snake25dList **slist, snake25d *snake) {
	snake25dList *newNode = (snake25dList *) malloc(sizeof(snake25dList));
	copySnake(snake, &newNode->snake);
	if (*slist) {
		newNode->next = *slist;
	} else {
		newNode->next = NULL;
	}
	*slist = newNode;
}

void removeSnake(snake25dList **slist) {
	snake25dList *node = *slist;
	if (node) {
		*slist = node->next;
		free(node);
	}
}

void destroySnakeList(snake25dList **slist) {
	while (*slist) {
		removeSnake(slist);
	}
}

int convertSnakeList2Array(snake25dList *slist, snake25d **sarray) {
	int n = 0;
	snake25dList *s = slist;
	while (s) {
		n++;
		s = s->next;
	}
	if (*sarray)
		free(*sarray);
	*sarray = NULL;
	if (n) {
		*sarray = (snake25d *) malloc(sizeof(snake25d) * n);
		n = 0;
		snake25dList *s = slist;
		while (s) {
			copySnake(&s->snake, &(*sarray)[n]);
			n++;
			s = s->next;
		}
	}
	return n;
}

void updateSnake(snake25d *s, double uv[][SN_N][2], const int w_1,
		const int h_1, double mov[][SN_N][2], double *maxMov) {
	int t, i;
	double px, py, dx, dy, m;
	for (t = 0; t < s->t; t++) {

		for (i = 0; i < SN_N; i++) {
			px = s->node[t][i][0];
			py = s->node[t][i][1];
			// Update snake
			s->node[t][i][0] += uv[t][i][0];
			s->node[t][i][1] += uv[t][i][1];
			// Viewpoint constraints & correction
			if (s->node[t][i][0] < 0)
				s->node[t][i][0] = 0;
			if (s->node[t][i][1] < 0)
				s->node[t][i][1] = 0;
			if (s->node[t][i][0] > w_1)
				s->node[t][i][0] = w_1;
			if (s->node[t][i][1] > h_1)
				s->node[t][i][1] = h_1;
			// Update movement
			dx = (mov[t][i][0] += s->node[t][i][0] - px);
			dy = (mov[t][i][1] += s->node[t][i][1] - py);
			m = dx * dx + dy * dy;
			if (*maxMov < m)
				*maxMov = m;
		}
	}
}

void resetMovement(double movement[][SN_N][2], double *maxMovement) {
	memset(movement, 0, sizeof(double) * SN25D_T * SN_N * 2);
	*maxMovement = 0;
}

void startSnake25d(curveEnergyGradientStack *grad_ces, snake25d *outputSnake,
		const double w_tens_xy, const double w_curv_xy, const double w_tens_z,
		const double w_curv_z, const double w_curve, const double w_inf,
		const double aspectRatio_z, snake25d *initialSnake, const double initx,
		const double inity, const int start_z, int end_z) {
	snake25d *snake = (snake25d *) malloc(sizeof(snake25d));
	double uv[SN25D_T][SN_N][2];
	double uv_corrected[SN25D_T][SN_N][2];
	double movement[SN25D_T][SN_N][2], maxMovement;

	int i, z;

	int w = grad_ces->grad_ces[0]->w - 1;
	int h = grad_ces->grad_ces[0]->h - 1;
	int w_1 = w - 1;
	int h_1 = h - 1;

	// Initiate snake
	if (initialSnake) {
		copySnake(initialSnake, snake);
	} else {
		initSnake25d(snake, initx, inity, start_z, end_z, SN25D_INITR, w, h,
				aspectRatio_z);
	}
	snake->status = SNAKE_GROWING;
	// Initialize iteration counter
	i = 0;
	// Reset cumulative snake movement
	resetMovement(movement, &maxMovement);
	// Main loop
	for (;;) {
		// Compute Snake Derivatives
		compSnakeDerivatives(snake);
		// Compute Internal Energy
		compIntEnergyGrad(snake, w_tens_xy, w_curv_xy, w_tens_z, w_curv_z);
		// Compute External Energy
		compExtEnergyGrad(snake, grad_ces, w_curve, w_inf);
		// Compute Update Vectors
		compUpdateVectors(snake, uv);
		// Equidistant correction
		equidistantCorrectionUV(snake, uv, uv_corrected);
		// Update snake
		updateSnake(snake, uv_corrected, w_1, h_1, movement, &maxMovement);
		// Compute area
		compSnakeArea(snake);
		// -----------------
		// Stopping criteria
		// -----------------
		// Max iter.
		if (i > SN25D_MAXITER) {
			snake->status = SNAKE_MAXITER;
			break;
		}
		// Min. area constraint
		for (z = 0; z < snake->t; z++) {
			if (snake->area[z] < SN25D_MINAREA) {
				snake->status = SNAKE_MINAREA;
				break;
			}
		}
		if (snake->status == SNAKE_MINAREA)
			break;
		// Max. area constraint
		for (z = 0; z < snake->t; z++) {
			if (snake->area[z] > SN25D_MAXAREA) {
				snake->status = SNAKE_MAXAREA;
				break;
			}
		}
		if (snake->status == SNAKE_MAXAREA)
			break;

		// Convergence criteria
		if (i % SN25D_SHORTTERM_CONV_ITER == SN25D_SHORTTERM_CONV_ITER - 1
				|| i % SN25D_LONGTERM_CONV_ITER
						== SN25D_LONGTERM_CONV_ITER - 1) {
			if (maxMovement < SN25D_SHORTTERM_CONV
					|| maxMovement < SN25D_LONGTERM_CONV) {
				snake->status = SNAKE_CONVERGED;
				break;
			}
			resetMovement(movement, &maxMovement);
		}
		// -----------------
		// Iterate
		i++;
	}
	// Return the output snake
	copySnake(snake, outputSnake);
	// Free mem.
	free(snake);
}

int inflationSnake25d(dataPacket *dp, const double initx, const double inity,
		const int start_z, int end_z, const double w_tens_xy,
		const double w_curv_xy, const double w_tens_z, const double w_curv_z,
		const double w_curve, const double w_inf_min, const double w_inf_max,
		const double w_inf_step, const double aspectRatio_z,
		snake25dList **outputSnakeList) {
	snake25d *snake = (snake25d *) malloc(sizeof(snake25d));
	snake25d *outputSnake = (snake25d *) malloc(sizeof(snake25d));
	double w_inf, sim;
	int flag = 0;
	int i = 0;

	for (w_inf = w_inf_min; w_inf <= w_inf_max; w_inf += w_inf_step) {
		startSnake25d(dp->grad_ces_lo, outputSnake, w_tens_xy, w_curv_xy,
				w_tens_z, w_curv_z, w_curve, w_inf, aspectRatio_z, NULL, initx,
				inity, start_z, end_z);
		if (outputSnake->status == SNAKE_CONVERGED) {
			// Check similarity
			sim = flag ? getDiceSimilaritySnake25d(snake, outputSnake) : 0.0;
			if (sim <= SN25D_INF_CONV) // If we have a different snake
			{
				// Validation
				validateSnake25d(dp, outputSnake);
				// Add to list
				pthread_mutex_lock(&addlock);
				addSnake(outputSnakeList, outputSnake);
				pthread_mutex_unlock(&addlock);
				// Report validity
				if (DT_REPORT) {
					reportValidation(dp, outputSnake);
				}
				i++;
				copySnake(outputSnake, snake);
				flag = 1; // Converged --save snake and restart with greater w_inf
			}
		}
		if (outputSnake->status == SNAKE_MAXAREA)
			break;	// Over inflation --stop
		if (outputSnake->status == SNAKE_MAXITER) {
			flag = 0;
		}				// Insufficient Inflation --restart with greater w_inf
		if (outputSnake->status == SNAKE_MINAREA) {
			flag = 0;
		}								// Vanished --restart with greater w_inf
	}

	// Free mem.
	free(snake);
	free(outputSnake);
	//
	return i;
}

void reportValidation(dataPacket *dp, snake25d *s) {
	for (int t = 0; t < s->t; t++) {
		isValidMitoSlice_25d(s, t, dp);
		cvWaitKey();
	}
}

void validateSnake25d(dataPacket *dp, snake25d *s) {
	int v, n_v = 0, n_iv = 0;
	for (int t = 0; t < s->t; t++) {
		v = isValidMitoSlice_25d(s, t, dp);
		if (v) {
			s->isValid[t] = SNAKE_VALID;
			n_v++;
		} else {
			s->isValid[t] = SNAKE_INVALID;
			n_iv++;
		}
	}
	s->validity = (double) n_v / (n_v + n_iv);
}

void validateSnakeList(dataPacket *dp, snake25dList *slist) {
	snake25dList *s = slist;
	while (s) {
		validateSnake25d(dp, &s->snake);
		s = s->next;
	}
}

int findInitialPointsWithinZRange(dataPacket *dp, CvPoint **pts, int start_z,
		int end_z) {
	int z, zci, zpi, ci, c;
	int i, j, x, y;
	double cx1, cx2, cy1, cy2, midx, midy, dx, dy, dm;
	double px, py;
	int n_avg, n_pts = 0, n_out;
	curve *tc;

	// Initialize memory
	int msize = 0;
	for (z = start_z; z <= end_z; z++)
		if (z >= dp->cs_lo->start_z && z <= dp->cs_lo->end_z)
			msize += dp->cs_lo->cstack[z - dp->cs_lo->start_z]->n;
	double (*rawPoints)[2] = (double (*)[2]) malloc(sizeof(double[2]) * msize);
	double (*avgPoints)[2] = (double (*)[2]) malloc(sizeof(double[2]) * msize);
	curve **assocCurve = (curve **) malloc(sizeof(curve *) * msize);
	int *neighbors = (int *) malloc(sizeof(int) * msize);
	char *visit = (char *) malloc(sizeof(char) * msize);
	memset(visit, 0, sizeof(char) * msize);
	double eps = SN25D_INITPTS_EPS;
	int minPts = SN25D_INITPTS_MIN;
	double distsq = eps * eps;
	double initloc = 2.0 * SN25D_INITR;

	// Determine the possible initial points
	for (z = start_z; z <= end_z; z++) {
		if (z >= dp->cs_lo->start_z && z <= dp->cs_lo->end_z
				&& z >= dp->cps_lo->start_z && z <= dp->cps_lo->end_z) {
			zci = z - dp->cs_lo->start_z;
			zpi = z - dp->cps_lo->start_z;
			ci = 0;
			for (c = 0; c < dp->cs_lo->cstack[zci]->n; c++) {
				// Store initial point
				if (dp->cs_lo->cstack[zci]->clist[c].h != 0) {
					cx1 = dp->cps_lo->cps[zpi]->cx[ci];
					cy1 = dp->cps_lo->cps[zpi]->cy[ci];
					cx2 = dp->cps_lo->cps[zpi]->cx[ci
							+ dp->cs_lo->cstack[zci]->clist[c].len - 1];
					cy2 = dp->cps_lo->cps[zpi]->cy[ci
							+ dp->cs_lo->cstack[zci]->clist[c].len - 1];
					midx = (cx1 + cx2) / 2.0;
					midy = (cy1 + cy2) / 2.0;
					dx = midx
							- dp->cps_lo->cps[zpi]->cx[ci
									+ dp->cs_lo->cstack[zci]->clist[c].len / 2];
					dy = midy
							- dp->cps_lo->cps[zpi]->cy[ci
									+ dp->cs_lo->cstack[zci]->clist[c].len / 2];
					dm = sqrt(dx * dx + dy * dy);
					dx *= initloc / dm;
					dy *= initloc / dm;
					rawPoints[n_pts][0] = dx
							+ dp->cps_lo->cps[zpi]->cx[ci
									+ dp->cs_lo->cstack[zci]->clist[c].len / 2];
					rawPoints[n_pts][1] = dy
							+ dp->cps_lo->cps[zpi]->cy[ci
									+ dp->cs_lo->cstack[zci]->clist[c].len / 2];
					assocCurve[n_pts] = &dp->cs_lo->cstack[zci]->clist[c];
					n_pts++;
				}

				ci += dp->cs_lo->cstack[zci]->clist[c].len;
			}
		}
	}

	// Density-based clustering
	for (i = 0; i < n_pts; i++) {
		c = 0;
		for (j = 0; j < n_pts; j++) {
			dx = rawPoints[i][0] - rawPoints[j][0];
			dy = rawPoints[i][1] - rawPoints[j][1];
			dm = dx * dx + dy * dy;
			if (dm <= distsq
					&& isInside(assocCurve[i], rawPoints[j][0], rawPoints[j][1])
					&& isInside(assocCurve[j], rawPoints[i][0],
							rawPoints[i][1])) {
				c++;
			}
		}
		neighbors[i] = c;
	}
	for (i = 0; i < n_pts - 1; i++) {
		for (j = i + 1; j < n_pts; j++) {
			if (neighbors[i] < neighbors[j]) {
				c = neighbors[i];
				neighbors[i] = neighbors[j];
				neighbors[j] = c;
				px = rawPoints[i][0];
				rawPoints[i][0] = rawPoints[j][0];
				rawPoints[j][0] = px;
				py = rawPoints[i][1];
				rawPoints[i][1] = rawPoints[j][1];
				rawPoints[j][1] = py;
				tc = assocCurve[i];
				assocCurve[i] = assocCurve[j];
				assocCurve[j] = tc;
			}
		}
	}
	n_avg = 0;
	for (i = 0; i < n_pts; i++) {
		if (!visit[i]) {
			c = 0;
			px = 0;
			py = 0;
			for (j = 0; j < n_pts; j++) {
				if (!visit[j]) {
					dx = rawPoints[i][0] - rawPoints[j][0];
					dy = rawPoints[i][1] - rawPoints[j][1];
					dm = dx * dx + dy * dy;
					if (dm <= distsq
							&& isInside(assocCurve[i], rawPoints[j][0],
									rawPoints[j][1])
							&& isInside(assocCurve[j], rawPoints[i][0],
									rawPoints[i][1])) {
						visit[j] = 1;
						px += rawPoints[j][0];
						py += rawPoints[j][1];
						c++;
					}
				}
			}
			if (c >= minPts) {
				avgPoints[n_avg][0] = px / c;
				avgPoints[n_avg][1] = py / c;
				n_avg++;
			}
		}
	}

	// Output results
	if (*pts)
		free(*pts);
	*pts = (CvPoint *) malloc(sizeof(CvPoint) * n_avg);
	n_out = 0;
	for (i = 0; i < n_avg; i++) {
		x = (int) avgPoints[i][0];
		y = (int) avgPoints[i][1];
		if (x > 0 && y > 0 && x < dp->ces_lo->ces[0]->w
				&& y < dp->ces_lo->ces[0]->h) {
			(*pts)[n_out].x = x;
			(*pts)[n_out].y = y;
			n_out++;
		}
	}

	// Free memory
	free(rawPoints);
	free(avgPoints);
	free(visit);
	free(assocCurve);

	return n_out;
}

struct t_param {
	dataPacket *dp;
	double ix, iy;
	int start_z, end_z;
	double ar_z;
	snake25dList **outputSnakeList;
	volatile int *returnValue;
};

void *threadPhase2(void *t_data) {
	struct t_param *p = (struct t_param *) t_data;
	int ns;
	ns = inflationSnake25d(p->dp, p->ix, p->iy, p->start_z, p->end_z,
	SN25D_W_TENSION,
	SN25D_W_CURVATURE,
	SN25D_W_ZTENSION,
	SN25D_W_ZCURVATURE,
	SN25D_W_ECURVE,
	SN25D_W_EINF_MIN,
	SN25D_W_EINF_MAX,
	SN25D_W_EINF_STEP, p->ar_z, p->outputSnakeList);

	pthread_mutex_lock(&tnslock);
	*p->returnValue += ns;
	pthread_mutex_unlock(&tnslock);
	sem_post(&sem2);
	pthread_exit(NULL);
}

int retrieveAllSnakes25d(dataPacket *dp, CvPoint *initPts, int n, int start_z,
		int end_z, snake25dList **outputSnakeList) {
	double ar_z = (int) LE_SSIZE_LO * TFACTOR / RESOLUTION;
	double ix, iy;
	int i;
	volatile int tns = 0;

	int t;
	int rc;
	int numThreads = n;
	pthread_t *threads = (pthread_t *) malloc(sizeof(pthread_t) * numThreads);
	pthread_attr_t attr;
	struct t_param *t_data = (struct t_param *) malloc(
			sizeof(struct t_param) * numThreads);
	sem_init(&sem2, 0, numCores);
	pthread_mutex_unlock(&addlock);
	pthread_mutex_unlock(&tnslock);
	void *status;
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

	// Report
	printf("Using %d initial points (slices: %04d - %04d)...\n", n, start_z,
			end_z);

	for (i = 0; i < n; i++) {
		sem_wait(&sem2);

		ix = (double) initPts[i].x;
		iy = (double) initPts[i].y;

		// Set thread params
		t_data[i].dp = dp;
		t_data[i].ix = ix;
		t_data[i].iy = iy;
		t_data[i].start_z = start_z;
		t_data[i].end_z = end_z;
		t_data[i].ar_z = ar_z;
		t_data[i].outputSnakeList = outputSnakeList;
		t_data[i].returnValue = &tns;
		//

		printf("\rProcessing: %d / %d (Total # of snakes : %d)", i + 1, n, tns);
		fflush(stdout);

		rc = pthread_create(&threads[i], &attr, threadPhase2,
				(void *) &t_data[i]);
		if (rc) {
			printf("Error: Unable to create thread!\n");
			exit(-1);
		}
	}

	pthread_attr_destroy(&attr);
	for (t = 0; t < numThreads; t++) {
		rc = pthread_join(threads[t], &status);
		if (rc) {
			printf("Error: Unable to join thread!\n");
			exit(-1);
		}
	}

	printf("\nTotal # of snakes : %d\n", tns);

	free(threads);
	free(t_data);

	return tns;
}

void saveSnakeArray(snake25d *sarray, int n, int tag_start_z, int tag_end_z) {
	FILE *f;
	char outputfn[256];
	char outputsrcfn[256];
	sprintf(outputsrcfn, "%s%s%d_%s.dat", DESTPATH, "snakes25d_",
			tag_end_z - tag_start_z + 1, FNAME);
	sprintf(outputfn, outputsrcfn, tag_start_z);
	printf("Saving: %s\n", outputfn);
	f = fopen(outputfn, "wb");

	fwrite(&n, sizeof(int), 1, f);
	for (int i = 0; i < n; i++) {
		fwrite(sarray[i].node, sizeof(double[SN_N][2]), sn25d_t, f);
		fwrite(sarray[i].d1_xy, sizeof(double[SN_N][2]), sn25d_t, f);
		fwrite(sarray[i].d2_xy, sizeof(double[SN_N][2]), sn25d_t, f);
		fwrite(sarray[i].d4_xy, sizeof(double[SN_N][2]), sn25d_t, f);
		fwrite(sarray[i].d1_z, sizeof(double[SN_N][2]), sn25d_t, f);
		fwrite(sarray[i].d2_z, sizeof(double[SN_N][2]), sn25d_t, f);
		fwrite(sarray[i].d4_z, sizeof(double[SN_N][2]), sn25d_t, f);
		fwrite(sarray[i].grad_Eint, sizeof(double[SN_N][2]), sn25d_t, f);
		fwrite(sarray[i].grad_Eext, sizeof(double[SN_N][2]), sn25d_t, f);
		fwrite(sarray[i].area, sizeof(double), sn25d_t, f);
		fwrite(&sarray[i].status, sizeof(char), 1, f);
		fwrite(sarray[i].isValid, sizeof(char), sn25d_t, f);
		fwrite(&sarray[i].validity, sizeof(double), 1, f);
		fwrite(&sarray[i].start_z, sizeof(int), 1, f);
		fwrite(&sarray[i].end_z, sizeof(int), 1, f);
		fwrite(&sarray[i].t, sizeof(int), 1, f);
		fwrite(&sarray[i].aspectRatio_z, sizeof(double), 1, f);
	}

	fclose(f);
}

int loadSnakeArray(int tag_start_z, int tag_end_z, snake25d **sarray) {
	int n = 0;
	FILE *f;
	size_t s;
	char outputfn[256];
	char outputsrcfn[256];
	sprintf(outputsrcfn, "%s%s%d_%s.dat", DESTPATH, "snakes25d_",
			tag_end_z - tag_start_z + 1, FNAME);
	sprintf(outputfn, outputsrcfn, tag_start_z);
	printf("Loading: %s\n", outputfn);
	f = fopen(outputfn, "rb");
	if (f) {
		s = fread(&n, sizeof(int), 1, f);
		if (*sarray)
			free(*sarray);
		*sarray = NULL;
		if (n > 0) {
			*sarray = (snake25d *) malloc(sizeof(snake25d) * n);
			//s = fread(*sarray, sizeof(snake25d), n, f);
			for (int i = 0; i < n; i++) {
				s = fread((*sarray)[i].node, sizeof(double[SN_N][2]), sn25d_t,
						f);
				s = fread((*sarray)[i].d1_xy, sizeof(double[SN_N][2]), sn25d_t,
						f);
				s = fread((*sarray)[i].d2_xy, sizeof(double[SN_N][2]), sn25d_t,
						f);
				s = fread((*sarray)[i].d4_xy, sizeof(double[SN_N][2]), sn25d_t,
						f);
				s = fread((*sarray)[i].d1_z, sizeof(double[SN_N][2]), sn25d_t,
						f);
				s = fread((*sarray)[i].d2_z, sizeof(double[SN_N][2]), sn25d_t,
						f);
				s = fread((*sarray)[i].d4_z, sizeof(double[SN_N][2]), sn25d_t,
						f);
				s = fread((*sarray)[i].grad_Eint, sizeof(double[SN_N][2]),
						sn25d_t, f);
				s = fread((*sarray)[i].grad_Eext, sizeof(double[SN_N][2]),
						sn25d_t, f);
				s = fread((*sarray)[i].area, sizeof(double), sn25d_t, f);
				s = fread(&((*sarray)[i].status), sizeof(char), 1, f);
				s = fread((*sarray)[i].isValid, sizeof(char), sn25d_t, f);
				s = fread(&((*sarray)[i].validity), sizeof(double), 1, f);
				s = fread(&((*sarray)[i].start_z), sizeof(int), 1, f);
				s = fread(&((*sarray)[i].end_z), sizeof(int), 1, f);
				s = fread(&((*sarray)[i].t), sizeof(int), 1, f);
				s = fread(&((*sarray)[i].aspectRatio_z), sizeof(double), 1, f);
			}
		}
		if (!s)
			printf("Error loading data!\n");
	} else
		printf("Error loading data!\n");
	fclose(f);
	return n;
}

void savePoints(CvPoint *pts, int n, int tag_start_z, int tag_end_z) {
	FILE *f;
	char outputfn[256];
	char outputsrcfn[256];
	sprintf(outputsrcfn, "%s%s%d_%s.dat", DESTPATH, "initPts_",
			tag_end_z - tag_start_z + 1, FNAME);
	sprintf(outputfn, outputsrcfn, tag_start_z);
	printf("Saving: %s\n", outputfn);
	f = fopen(outputfn, "wb");

	fwrite(&n, sizeof(int), 1, f);
	fwrite(pts, sizeof(CvPoint), n, f);

	fclose(f);
}

int loadPoints(int tag_start_z, int tag_end_z, CvPoint **pts) {
	int n = 0;
	FILE *f;
	size_t s;
	char outputfn[256];
	char outputsrcfn[256];
	sprintf(outputsrcfn, "%s%s%d_%s.dat", DESTPATH, "initPts_",
			tag_end_z - tag_start_z + 1, FNAME);
	sprintf(outputfn, outputsrcfn, tag_start_z);
	printf("Loading: %s\n", outputfn);
	f = fopen(outputfn, "rb");
	if (f) {
		s = fread(&n, sizeof(int), 1, f);
		if (*pts)
			free(*pts);
		*pts = NULL;
		if (n > 0) {
			*pts = (CvPoint *) malloc(sizeof(CvPoint) * n);
			s = fread(*pts, sizeof(CvPoint), n, f);
		}
		if (!s)
			printf("Error loading data!\n");
	} else
		printf("Error loading data!\n");
	fclose(f);
	return n;
}

int filterSnakeArrayByValidity(snake25d *in, int n, double th, snake25d **out) {
	int i, j, n_out = 0;
	for (i = 0; i < n; i++)
		if (in[i].validity >= th)
			n_out++;
	if (*out)
		free(*out);
	*out = NULL;
	if (n_out) {
		*out = (snake25d *) malloc(sizeof(snake25d) * n_out);
		j = 0;
		for (i = 0; i < n; i++) {
			if (in[i].validity >= th) {
				(*out)[j] = in[i];
				j++;
			}
		}
	}
	return n_out;
}

void filterSnake25dByValidity(snake25dList *in, snake25dList **out, double th) {
	if (*out)
		destroySnakeList(out);
	*out = NULL;
	snake25dList *s = in;
	while (s) {
		if (s->snake.validity >= th) {
			addSnake(out, &s->snake);
		}
		s = s->next;
	}
}

int disintegrateSnake25dTo2d(snake25d *in, snake25d *outArray) {
	int i;
	for (i = 0; i < in->t; i++) {
		outArray[i].area[0] = in->area[i];
		outArray[i].aspectRatio_z = in->aspectRatio_z;
		memcpy(outArray[i].d1_xy[0], in->d1_xy[i], sizeof(double[SN_N][2]));
		memset(outArray[i].d1_z[0], 0, sizeof(double[SN_N][2]));
		memcpy(outArray[i].d2_xy[0], in->d2_xy[i], sizeof(double[SN_N][2]));
		memset(outArray[i].d2_z[0], 0, sizeof(double[SN_N][2]));
		memcpy(outArray[i].d4_xy[0], in->d4_xy[i], sizeof(double[SN_N][2]));
		memset(outArray[i].d4_z[0], 0, sizeof(double[SN_N][2]));
		outArray[i].start_z = in->start_z + i;
		outArray[i].end_z = outArray[i].start_z;
		outArray[i].t = 1;
		memcpy(outArray[i].grad_Eext[0], in->grad_Eext[i],
				sizeof(double[SN_N][2]));
		memcpy(outArray[i].grad_Eint[0], in->grad_Eint[i],
				sizeof(double[SN_N][2]));
		outArray[i].isValid[0] = in->isValid[i];
		memcpy(outArray[i].node[0], in->node[i], sizeof(double[SN_N][2]));
		outArray[i].status = in->status;
		outArray[i].validity = in->isValid[i] ? 1.0 : 0.0;
	}
	return in->t;
}

int disintegrateArrayOfSnake25dTo2d(snake25d *sarray, int n,
		snake25d **outArray) {
	int n_snk = 0, i, t = 0;
	for (i = 0; i < n; i++)
		n_snk += sarray[i].t;

	if (*outArray)
		free(*outArray);
	*outArray = NULL;
	if (n_snk) {
		*outArray = (snake25d *) malloc(sizeof(snake25d) * n_snk);

		for (i = 0; i < n; i++) {
			t += disintegrateSnake25dTo2d(&sarray[i], &(*outArray)[t]);
		}
	}
	return n_snk;
}

int filterArrayOfSnake25dBySimilarity(snake25d *sarray, int n,
		snake25d **outArray, double th) {
	int i, j;
	int z1, z2;
	double ds;
	char *filter;
	int n_out = 0;

	if (n) {
		filter = (char *) malloc(sizeof(char) * n);
		memset(filter, 0, sizeof(char) * n);

		n_out = n;
		for (i = 0; i < n - 1; i++) {
			for (j = i + 1; j < n; j++) {
				if (filter[j])
					continue;
				z1 = (sarray[i].start_z > sarray[j].start_z) ?
						sarray[i].start_z : sarray[j].start_z;
				z2 = (sarray[i].end_z < sarray[j].end_z) ?
						sarray[i].end_z : sarray[j].end_z;
				if (z1 <= z2) {
					ds = getDiceSimilaritySnake25d(&sarray[i], &sarray[j]);
					if (ds >= th) {
						if (sarray[i].validity > sarray[j].validity) {
							filter[j] = 1;
							n_out--;
						} else {
							filter[i] = 1;
							n_out--;
							break;
						}
					}
				}
			}
		}

		if (*outArray)
			free(*outArray);
		*outArray = NULL;
		if (n_out) {
			*outArray = (snake25d *) malloc(sizeof(snake25d) * n_out);

			j = 0;
			for (i = 0; i < n; i++) {
				if (!filter[i]) {
					(*outArray)[j] = sarray[i];
					j++;
				}
			}
		}

		free(filter);
	}

	return n_out;
}

void getEnergyGradient(const double *ce, const int cw, const int ch,
		double (*grad_ce)[2]) {
	int x, y;
	int xp, yp;
	double e_xp, e_yp, e;
	for (y = 0; y < ch; y++) {
		yp = y + 1;
		for (x = 0; x < cw; x++) {
			e = ce[y * cw + x];
			xp = x + 1;
			e_xp = e_yp = 0;
			if (xp < cw)
				e_xp = ce[y * cw + xp];
			if (yp < ch)
				e_yp = ce[yp * cw + x];
			grad_ce[y * cw + x][0] = e_xp - e;
			grad_ce[y * cw + x][1] = e_yp - e;
		}
	}
}

void getCurveEnergyImage(curve *clist, const int n, const int w, const int h,
		double **ce) {
	int i, j;
	curve *c;
	double *cx, *cy;
	IplImage *ground;
	double max1, max2, min1, min2;
	int msize;
	int xi, yi;

	msize = w * h;

	// Dealloc mem if needed
	if (*ce)
		free(*ce);

	// Alloc mem
	ground = cvCreateImage(cvSize(w, h), IPL_DEPTH_32F, 1);
	cvSetZero(ground);

	// Compute ground curve energy
	for (i = 0; i < n; i++) {
		c = &clist[i];

		cx = (double *) malloc(sizeof(double) * c->len);
		cy = (double *) malloc(sizeof(double) * c->len);

		getCurvePoints(c, 1.0, cx, cy);

		for (j = 0; j < c->len; j++) {
			xi = cvRound(cx[j]);
			yi = cvRound(cy[j]);
			if (xi >= 0 && yi >= 0 && xi < w && yi < h)
				cvSetReal2D(ground, yi, xi,
						cvGetReal2D(ground, yi, xi) + c->avgscore);
		}

		free(cx);
		free(cy);
	}

	// Smoothing
	cvMinMaxLoc(ground, &min1, &max1);
	cvSmooth(ground, ground, CV_GAUSSIAN, 0, 0, SN_GAUSSIAN);
	cvMinMaxLoc(ground, &min2, &max2);
	cvConvertScale(ground, ground, max1 / max2);

	// Alloc mem
	*ce = (double *) malloc(sizeof(double) * msize);
	memset(*ce, 0, sizeof(double) * msize);

	// Copy to output array
	for (i = 0, yi = 0; yi < h; yi++) {
		for (xi = 0; xi < w; xi++, i++) {
			(*ce)[i] = cvGetReal2D(ground, yi, xi);
		}
	}

	// Free mem.
	cvReleaseImage(&ground);
}
