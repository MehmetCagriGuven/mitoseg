/*
 * datasets.h
 *
 *  Created on: Jan 13, 2017
 *      Author: fst
 */

#ifndef DATASETS_H_
#define DATASETS_H_

#include "opencv/cv.h"
#include "opencv/highgui.h"

extern char FNAME[256];
extern char DESTPATH[256];
extern char SOURCEPATH[256];
extern int SLICE_START;
extern int SLICE_END;
extern int ROI_X;
extern int ROI_Y;
extern int ROI_W;
extern int ROI_H;
extern double RESOLUTION;

int loadSlice(int s, IplImage **im);
void saveImage(int s, char const *tag, char const *tag2, IplImage *im);
void saveImageDat(int s, char const *tag, IplImage *im);
void loadImageDat(int s, char const *tag, IplImage **im);
void getImageROI(IplImage *input, IplImage **output);

#endif /* DATASETS_H_ */
