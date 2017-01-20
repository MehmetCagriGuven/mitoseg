/*
 * datasets.cpp
 *
 *  Created on: Jan 13, 2017
 *      Author: fst
 */

#include "datasets.h"

#include <stdio.h>

char FNAME[256];
char DESTPATH[256];
char SOURCEPATH[256];
int SLICE_START;
int SLICE_END;
int ROI_X;
int ROI_Y;
int ROI_W;
int ROI_H;
double RESOLUTION;

int loadSlice(int s, IplImage **im) {
	char inputfn[256];
	char inputsrcfn[256];
	sprintf(inputsrcfn, "%s%s", SOURCEPATH, FNAME);
	sprintf(inputfn, inputsrcfn, s);
	printf("Loading: %s\n", inputfn);
	if (*im)
		cvReleaseImage(im);
	*im = cvLoadImage(inputfn);
	if (!(*im)) {
		printf("Error loading image!\n");
		return 1;
	}
	return 0;
}

void saveImage(int s, char const *tag, char const *tag2, IplImage *im) {
	char outputfn[256];
	char outputsrcfn[256];
	sprintf(outputsrcfn, "%s%s%s%s", DESTPATH, tag, tag2, FNAME);
	sprintf(outputfn, outputsrcfn, s);
	if (im) {
		printf("Saving: %s\n", outputfn);
		cvSaveImage(outputfn, im);
	}
}

void saveImageDat(int s, char const *tag, IplImage *im) {
	char outputfn[256];
	char outputsrcfn[256];
	sprintf(outputsrcfn, "%s%s%s.dat", DESTPATH, tag, FNAME);
	sprintf(outputfn, outputsrcfn, s);
	if (im) {
		printf("Saving: %s\n", outputfn);
		cvSave(outputfn, im);
	}
}

void loadImageDat(int s, char const *tag, IplImage **im) {
	char outputfn[256];
	char outputsrcfn[256];
	sprintf(outputsrcfn, "%s%s%s.dat", DESTPATH, tag, FNAME);
	sprintf(outputfn, outputsrcfn, s);
	if (*im) {
		cvReleaseImage(im);
	}
	printf("Loading: %s\n", outputfn);
	*im = (IplImage *) cvLoad(outputfn);
}

CvScalar operator*(const double &c, const CvScalar &a) {
	CvScalar t;
	t.val[0] = a.val[0] * c;
	t.val[1] = a.val[1] * c;
	t.val[2] = a.val[2] * c;
	t.val[3] = a.val[3] * c;
	return t;
}

CvScalar operator-(const CvScalar &a, const CvScalar &b) {
	CvScalar t;
	t.val[0] = a.val[0] - b.val[0];
	t.val[1] = a.val[1] - b.val[1];
	t.val[2] = a.val[2] - b.val[2];
	t.val[3] = a.val[3] - b.val[3];
	return t;
}

double operator~(const CvScalar &a) {
	double d;
	d = a.val[0] * a.val[0];
	d += a.val[1] * a.val[1];
	d += a.val[2] * a.val[2];
	d += a.val[3] * a.val[3];
	return d;
}

double getDev(IplImage *input) {
	double d = 0;
	for (int y = 1; y < input->height - 1; y++) {
		for (int x = 1; x < input->width - 1; x++) {
			d += ~(0.5 * cvGet2D(input, y, x) - 0.125 * cvGet2D(input, y - 1, x)
					- 0.125 * cvGet2D(input, y + 1, x)
					- 0.125 * cvGet2D(input, y, x - 1)
					- 0.125 * cvGet2D(input, y, x + 1));
		}
	}
	return d / (input->height * input->width);
}

bool checkRow(IplImage *input, const int y, const double &dev) {
	double d = 0;
	int yp, ym;
	yp = ym = y;
	if (y > 0)
		ym--;
	if (y < input->height - 1)
		yp++;
	for (int x = 0; x < input->width; x++)
		d += ~(0.5 * cvGet2D(input, y, x) - 0.25 * cvGet2D(input, yp, x)
				- 0.25 * cvGet2D(input, ym, x));
	d /= input->width;
	return 6.0 * d < dev;
}

bool checkColumn(IplImage *input, const int x, const double &dev) {
	double d = 0;
	int xp, xm;
	xp = xm = x;
	if (x > 0)
		xm--;
	if (x < input->width - 1)
		xp++;
	for (int y = 0; y < input->height; y++)
		d += ~(0.5 * cvGet2D(input, y, x) - 0.25 * cvGet2D(input, y, xp)
				- 0.25 * cvGet2D(input, y, xm));
	d /= input->width;
	return 6.0 * d < dev;
}

void autoROI(IplImage *input) {
	double dev = getDev(input);
	ROI_X = ROI_Y = 0;
	ROI_W = input->width;
	ROI_H = input->height;
	while (ROI_H > 1 && checkRow(input, ROI_Y, dev)) {
		ROI_Y++;
		ROI_H--;
	}
	while (ROI_H > 1 && checkRow(input, ROI_Y + ROI_H - 1, dev)) {
		ROI_H--;
	}
	while (ROI_W > 1 && checkColumn(input, ROI_X, dev)) {
		ROI_X++;
		ROI_W--;
	}
	while (ROI_W > 1 && checkColumn(input, ROI_X + ROI_W - 1, dev)) {
		ROI_W--;
	}
	ROI_X += 5;
	ROI_Y += 5;
	ROI_W -= 10;
	ROI_H -= 10;
}

void getImageROI(IplImage *input, IplImage **output) {
	if (*output)
		cvReleaseImage(output);
	if (ROI_X < 0) {
		autoROI(input);
		printf("Using ROI: %d %d %d %d\n", ROI_X, ROI_Y, ROI_W, ROI_H);
	}
	*output = cvCreateImage(cvSize(ROI_W, ROI_H), input->depth,
			input->nChannels);
	cvSetImageROI(input, cvRect(ROI_X, ROI_Y, ROI_W, ROI_H));
	cvCopy(input, *output);
	cvResetImageROI(input);
}
