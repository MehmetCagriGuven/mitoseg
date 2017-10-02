/*
 * poly.h
 *
 *  Created on: Jan 15, 2017
 *      Author: fst
 */

#ifndef POLY_H_
#define POLY_H_

#include "snake25d.h"

struct poly2d {
	CvPoint *vertex;
	int n;
	poly2d *next;
};

struct poly25d {
	poly2d *slice;
	int t;
	int start_z, end_z;
};

void initBlankPoly2d(int n, poly2d *p);
void initPoly2d(double (*v)[2], int n, poly2d *p);
void deinitPoly2d(poly2d *p);
void initBlankPoly25d(int t, int start_z, poly25d *p);
void deinitPoly25d(poly25d *p);
void initPoly25d(snake25d *s, poly25d *p);
void initPoly25dArray(snake25d *sarray, int n, poly25d *parray);
int poly2dArea(poly2d *p, int *r = NULL);
void correctPoly2d(poly2d *p);
void deinitPoly25dArray(poly25d *parray, int n);
void combinePoly2d(poly2d *p1, poly2d *p2, int intersect = 0,
		poly2d *out = NULL, int *a1 = NULL, int *a2 = NULL, int *a12 = NULL);
void combinePoly25d(poly25d *p1, poly25d *p2, int intersect = 0, poly25d *out =
NULL, int *a1 = NULL, int *a2 = NULL, int *a12 = NULL);
void safecopyPoly2d(poly2d *source, poly2d *dest);
void safecopyPoly25d(poly25d *source, poly25d *dest);
int mergeArrayOfPoly25d(poly25d *parray, int n, double th, poly25d *outArray);
int convertValidSnakesToPolyArray(snake25d *sarray, int n, double th,
		poly25d **parray);
void savePolyArrayAsPLY(poly25d *p, int n);
int getTriangles(poly2d *p1, poly2d *p2, int (**faces)[3]);

#endif /* POLY_H_ */
