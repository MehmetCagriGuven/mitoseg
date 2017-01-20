/*
 * visualizer.h
 *
 *  Created on: Jan 15, 2017
 *      Author: fst
 */

#ifndef VISUALIZER_H_
#define VISUALIZER_H_

#include "opencv/cv.h"
#include "opencv/highgui.h"

#include "datasets.h"
#include "snake25d.h"
#include "poly.h"

class Visualizer {

public:
	IplImage *outputImage;
	//
	Visualizer();
	~Visualizer();
	void setStartSlice(int s);
	int getStartSlice();
	void setNumSlices(int n);
	int getNumSlices();
	void setFileNameTag(const char *path, const char *tag1, const char *tag2);
	void update();
	void flushBuffer();
	void getSliceCoor(CvPoint p, int *s, CvPoint *c);
	void setSnakeBuf(snake25d *buf, int n = 1, double scale = 1.0);
	void setMarkerBuf(CvPoint *buf, int n = 0, int z = -1, double scale = 1.0);
	void setPolyBuf(poly25d *buf, int n = 1, double scale = 1.0);
	CvSize getFrameSize();
	CvSize getGridSize();
	CvSize getOutputSize();

private:
	int startSlice;
	int numSlices;
	int w, h, gx, gy, gw, gh, width, height;
	char fname_path[256], fname_tag1[256], fname_tag2[256];
	CvSize sliceSize;
	//
	IplImage *buf[MAXBUFFER];
	int bufIndex[MAXBUFFER];
	int nextFree;
	//
	snake25d *snakeBuf;
	int n_snakeBuf;
	double snake_scale;
	//
	CvPoint *markerBuf;
	int n_markerBuf;
	double marker_scale;
	int marker_z;
	//
	poly25d *polyBuf;
	int n_polyBuf;
	double poly_scale;
	//
	int getBuffer(int s, IplImage **im);
	void insertBuffer(int s, IplImage *im);
	int loadSliceImage(int s, IplImage **im);
	IplImage* getSlice(int slice);
	IplImage* getFrame(int slice, CvSize size);
	void getNATemplate(IplImage **buffer);
	void drawFrame(int slice, IplImage *buffer);
	void drawSnake(IplImage *buffer, snake25d *snake, int snakeSlice);
	void drawMarker(IplImage *buffer);
	void drawPoly(IplImage *buffer, poly25d *poly, int polySlice);
};

#endif /* VISUALIZER_H_ */
