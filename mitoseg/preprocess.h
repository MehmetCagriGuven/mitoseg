/*
 * preprocess.h
 *
 *  Created on: Jan 13, 2017
 *      Author: fst
 */

#ifndef PREPROCESS_H_
#define PREPROCESS_H_

#include "opencv/cv.h"
#include "opencv/highgui.h"

#include <stdio.h>

#include "datasets.h"
#include "settings.h"
#include "definitions.h"

// Ridge Structure
typedef struct {
	float angle;
	float dir;
	int dir4;
	int dir8;
	float mag;
	float vec[2];
	float val[2];
	char mask;
} ridge;
//////////

////////////
void getResampledImage(IplImage *input, IplImage **output);
void autoLevels(IplImage *input, IplImage *output);
void getRidges(IplImage *input, ridge *r, IplImage *output,
		IplImage *output_dir, IplImage *output_dir4, IplImage *output_dir8,
		IplImage *mask);
void visualizeRidgesDir(int opt, ridge *r, IplImage *visual_dir);
void normalizeRidges(ridge *r, IplImage *ridges, IplImage *mask);
void fastBilateralFilter(IplImage *input, IplImage *output, int spatialWindow,
		float graySigma, int order);
void saveRidgeStructure(int s, char const *tag, ridge *r, int numRecords);
int loadRidgeStructure(int s, char const *tag, ridge **r);
/////////////

#endif /* PREPROCESS_H_ */
