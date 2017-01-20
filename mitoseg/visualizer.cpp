/*
 * visualizer.cpp
 *
 *  Created on: Jan 15, 2017
 *      Author: fst
 */

#include "visualizer.h"

#include <stdio.h>
#include <string.h>

void visualizeSnake(IplImage *visual, const double snake[][2],
		const double scale, CvScalar cl = cvScalar(0, 255, 255),
		int markerSize = 3, CvScalar cl2 = cvScalar(255, 255, 255)) {
	int i, j, x, y, xn, yn;
	for (i = 0; i < SN_N; i++) {
		j = (i + 1) % SN_N;
		x = cvRound(scale * snake[i][0]);
		y = cvRound(scale * snake[i][1]);
		xn = cvRound(scale * snake[j][0]);
		yn = cvRound(scale * snake[j][1]);
		cvLine(visual, cvPoint(x, y), cvPoint(xn, yn), cl2);
		cvCircle(visual, cvPoint(x, y), markerSize, cl, CV_FILLED);
	}
}

Visualizer::Visualizer() {
	outputImage = NULL;
	//
	startSlice = 0;
	numSlices = 1;
	//
	gx = gy = 1;
	gw = gh = w = h = width = height = 0;
	//
	for (int i = 0; i < MAXBUFFER; i++) {
		buf[i] = NULL;
		bufIndex[i] = -1;
	}
	nextFree = 0;
	//
	setFileNameTag("", "", "");
	//
	snakeBuf = NULL;
	n_snakeBuf = 0;
	snake_scale = 1.0;
	//
	markerBuf = NULL;
	n_markerBuf = 0;
	marker_z = -1;
	marker_scale = 1.0;
	//
	polyBuf = NULL;
	n_polyBuf = 0;
	poly_scale = 1.0;
}

Visualizer::~Visualizer() {
	if (outputImage)
		cvReleaseImage(&outputImage);
	flushBuffer();
}

void Visualizer::setSnakeBuf(snake25d *buf, int n, double scale) {
	snakeBuf = buf;
	n_snakeBuf = n;
	snake_scale = scale;
}

void Visualizer::setMarkerBuf(CvPoint *buf, int n, int z, double scale) {
	markerBuf = buf;
	n_markerBuf = n;
	marker_z = z;
	marker_scale = scale;
}

void Visualizer::setPolyBuf(poly25d *buf, int n, double scale) {
	polyBuf = buf;
	n_polyBuf = n;
	poly_scale = scale;
}

void Visualizer::setStartSlice(int s) {
	if (s >= SLICE_START && s <= SLICE_END)
		startSlice = s;
}

int Visualizer::getStartSlice() {
	return startSlice;
}

void Visualizer::setNumSlices(int n) {
	if (n >= 1 && n <= MAXIMAGES)
		numSlices = n;
}

int Visualizer::getNumSlices() {
	return numSlices;
}

CvSize Visualizer::getFrameSize() {
	return cvSize(width, height);
}

CvSize Visualizer::getGridSize() {
	return cvSize(gw, gh);
}

CvSize Visualizer::getOutputSize() {
	return cvSize(gx, gy);
}

void Visualizer::setFileNameTag(const char *path, const char *tag1,
		const char *tag2) {
	strcpy(fname_tag1, tag1);
	strcpy(fname_tag2, tag2);
	strcpy(fname_path, path);
	flushBuffer();
	if (fname_tag1[0] == 0 && fname_tag2[0] == 0 && fname_path[0] == 0) {
		sliceSize = cvSize(ROI_W, ROI_H);
	} else {
		IplImage *buffer = NULL;
		for (int s = SLICE_START; s <= SLICE_END; s++) {
			buffer = getSlice(s);
			if (buffer) {
				sliceSize = cvSize(buffer->width, buffer->height);
				cvReleaseImage(&buffer);
				return;
			}
		}

		if (!buffer)
			sliceSize = cvSize(50, 50);
	}
}

int Visualizer::getBuffer(int s, IplImage **im) {
	for (int i = 0; i < MAXBUFFER; i++) {
		if (bufIndex[i] == s) {
			*im = buf[i];
			return 0;
		}
	}
	*im = NULL;
	return -1;
}

void Visualizer::insertBuffer(int s, IplImage *im) {
	int i, f = 0, ind = nextFree;
	for (i = 0; i < MAXBUFFER; i++) {
		if (bufIndex[i] == s) {
			ind = i;
			f = 1;
			break;
		}
	}

	if (buf[ind])
		cvReleaseImage(&buf[ind]);
	buf[ind] = cvCloneImage(im);
	bufIndex[ind] = s;

	if (!f)
		nextFree = (nextFree + 1) % MAXBUFFER;
}

int Visualizer::loadSliceImage(int s, IplImage **im) {
	char fname[256], destpath[256], fn[256];
	IplImage *temp = NULL;

	if (fname_tag1[0] == 0 && fname_tag2[0] == 0 && fname_path[0] == 0) {
		int r = loadSlice(s, &temp);
		if (!r) {
			getImageROI(temp, im);
			cvReleaseImage(&temp);
		}
		return r;
	} else {
		if (fname_path[0] != 0)
			strcpy(destpath, fname_path);
		else
			strcpy(destpath, DESTPATH);
		sprintf(fn, "%s%s%s%s", destpath, fname_tag1, fname_tag2, FNAME);
		sprintf(fname, fn, s);
		printf("Loading: %s\n", fname);
		if (*im)
			cvReleaseImage(im);
		*im = cvLoadImage(fname);
		if ((*im) == NULL) {
			printf("Error loading image!\n");
			return 1;
		}
		if ((*im)->nChannels == 1) {
			temp = (*im);
			*im = cvCreateImage(cvSize(temp->width, temp->height), temp->depth,
					3);
			cvMerge(temp, temp, temp, NULL, *im);
			cvReleaseImage(&temp);
		}
		return 0;
	}
}

IplImage* Visualizer::getSlice(int slice) {
	IplImage *buffer = NULL;
	if (slice >= SLICE_START && slice <= SLICE_END) {
		loadSliceImage(slice, &buffer);
	}
	return buffer;
}

IplImage* Visualizer::getFrame(int slice, CvSize size) {
	IplImage *temp = NULL, *buffer = NULL;
	if (!getBuffer(slice, &buffer)) {
		if (buffer->width == size.width && buffer->height == size.height)
			return cvCloneImage(buffer);
	}

	buffer = getSlice(slice);
	if (buffer) {
		temp = cvCreateImage(size, buffer->depth, buffer->nChannels);
		cvResize(buffer, temp, CV_INTER_AREA);
		cvReleaseImage(&buffer);
		insertBuffer(slice, temp);
		return temp;
	}

	getNATemplate(&buffer);
	temp = cvCreateImage(size, buffer->depth, buffer->nChannels);
	cvResize(buffer, temp, CV_INTER_NN);
	cvReleaseImage(&buffer);
	return temp;
}

void Visualizer::getNATemplate(IplImage **buffer) {
	CvFont font;
	CvSize text_size;
	int baseline;

	if (*buffer)
		cvReleaseImage(buffer);
	*buffer = cvCreateImage(cvSize(50, 50), 8, 3);
	cvSetZero(*buffer);
	cvInitFont(&font, CV_FONT_HERSHEY_SIMPLEX, 0.3, 0.3);
	cvGetTextSize("N/A", &font, &text_size, &baseline);
	cvPutText(*buffer, "N/A",
			cvPoint(25 - text_size.width / 2, 25 - text_size.height / 2), &font,
			cvScalar(0, 255, 255));
}

void Visualizer::drawSnake(IplImage *buffer, snake25d *snake, int snakeSlice) {
	CvScalar cl;
	int marker;
	switch (snake->status) {
	case SNAKE_GROWING:
		cl = cvScalar(0, 255, 255);
		break;
	case SNAKE_CONVERGED:
		if (snake->isValid[snakeSlice] == SNAKE_VALID)
			cl = cvScalar(0, 255, 0);
		else
			cl = cvScalar(0, 0, 255);
		break;
	case SNAKE_MAXAREA:
	case SNAKE_MINAREA:
	case SNAKE_MAXITER:
		cl = cvScalar(0, 0, 255);
		break;
	default:
		cl = cvScalar(255, 255, 255);
	}
	marker = MAX(cvRound(5.0 * gw / width), 1);
	visualizeSnake(buffer, snake->node[snakeSlice], snake_scale * gw / width,
			cl, marker);
}

void Visualizer::drawMarker(IplImage *buffer) {
	int i;
	int x, y;
	int s = MAX(cvRound(10.0 * gw / width), 1);
	int t = MAX(s / 2, 1);
	for (i = 0; i < n_markerBuf; i++) {
		x = cvRound(markerBuf[i].x * marker_scale * gw / width);
		y = cvRound(markerBuf[i].y * marker_scale * gw / width);
		cvLine(buffer, cvPoint(x - s, y), cvPoint(x + s, y),
				cvScalar(0, 255, 255), t);
		cvLine(buffer, cvPoint(x, y - s), cvPoint(x, y + s),
				cvScalar(0, 255, 255), t);
	}
}

void Visualizer::drawPoly(IplImage *buffer, poly25d *poly, int polySlice) {
	int t = MAX(cvRound(5.0 * gw / width), 1);
	int x1, y1, x2, y2, i;
	poly2d *seq = &poly->slice[polySlice];
	while (seq) {
		for (i = 0; i < seq->n; i++) {
			x1 = cvRound(seq->vertex[i].x * poly_scale * gw / width);
			y1 = cvRound(seq->vertex[i].y * poly_scale * gw / width);
			x2 = cvRound(
					seq->vertex[(i + 1) % seq->n].x * poly_scale * gw / width);
			y2 = cvRound(
					seq->vertex[(i + 1) % seq->n].y * poly_scale * gw / width);
			cvLine(buffer, cvPoint(x1, y1), cvPoint(x2, y2),
					cvScalar(255, 0, 0), t);
		}
		seq = seq->next;
	}
}

void Visualizer::drawFrame(int slice, IplImage *buffer) {
	int n;
	// Draw snakes
	if (snakeBuf) {
		int ss;
		for (n = 0; n < n_snakeBuf; n++) {
			if (slice >= snakeBuf[n].start_z && slice <= snakeBuf[n].end_z) {
				ss = slice - snakeBuf[n].start_z;
				drawSnake(buffer, &snakeBuf[n], ss);
			}
		}
	}
	// Draw markers
	if (markerBuf && slice == marker_z) {
		drawMarker(buffer);
	}
	// Draw poly
	if (polyBuf) {
		int ss;
		for (n = 0; n < n_polyBuf; n++) {
			if (slice >= polyBuf[n].start_z && slice <= polyBuf[n].end_z) {
				ss = slice - polyBuf[n].start_z;
				drawPoly(buffer, &polyBuf[n], ss);
			}
		}
	}
}

void Visualizer::flushBuffer() {
	int i;
	for (i = 0; i < MAXBUFFER; i++) {
		if (buf[i])
			cvReleaseImage(&buf[i]);
		buf[i] = NULL;
		bufIndex[i] = -1;
	}
	nextFree = 0;
}

void Visualizer::update() {
	IplImage *buffer = NULL;
	CvFont font;
	CvSize text_size;
	int baseline;

	width = sliceSize.width;
	height = sliceSize.height;

	gx = gy = 1;
	double r = (double) WINDOW_WIDTH / WINDOW_HEIGHT;
	double ri = (double) width / height;
	while (gx * gy < numSlices) {
		double rgx = ri * (gx + 1) / gy;
		double rgy = ri * gx / (gy + 1);
		if (fabs(rgx - r) < fabs(rgy - r))
			gx++;
		else
			gy++;
	}

	w = gx * width;
	h = gy * height;
	gw = width;
	gh = height;
	if (w > WINDOW_WIDTH || h > WINDOW_HEIGHT) {
		int _h = h * WINDOW_WIDTH / w;
		if (_h <= WINDOW_HEIGHT) {
			gw = gw * WINDOW_WIDTH / w;
			gh = gh * WINDOW_WIDTH / w;
		} else {
			gw = gw * WINDOW_HEIGHT / h;
			gh = gh * WINDOW_HEIGHT / h;
		}
	}
	w = gx * gw;
	h = gy * gh;

	if (outputImage)
		cvReleaseImage(&outputImage);
	outputImage = cvCreateImage(cvSize(w, h), 8, 3);
	cvSetZero(outputImage);

	cvInitFont(&font, CV_FONT_HERSHEY_SIMPLEX, 0.3, 0.3);
	char fn[10];

	int x, y, px, py, i, s;
	for (i = 0, y = 0, py = 0, s = startSlice;
			y < gy && i < numSlices && s <= SLICE_END; y++, py += gh) {
		for (x = 0, px = 0; x < gx && i < numSlices && s <= SLICE_END;
				x++, px += gw, i++, s++) {
			// Frame
			buffer = getFrame(s, cvSize(gw, gh));
			drawFrame(s, buffer);
			cvSetImageROI(outputImage, cvRect(px, py, gw, gh));
			cvCopy(buffer, outputImage);
			cvResetImageROI(outputImage);
			cvReleaseImage(&buffer);
			// Frame #
			sprintf(fn, "%d", i + startSlice);
			cvGetTextSize(fn, &font, &text_size, &baseline);
			cvPutText(outputImage, fn,
					cvPoint(px + 2, py + text_size.height + 2), &font,
					cvScalar(0, 255, 255));
		}
	}
	// Post drawings
	for (y = 1, py = gh; y < gy; y++, py += gh) {
		cvLine(outputImage, cvPoint(0, py), cvPoint(w - 1, py),
				cvScalar(0, 0, 255));
	}
	for (x = 1, px = gw; x < gx; x++, px += gw) {
		cvLine(outputImage, cvPoint(px, 0), cvPoint(px, h - 1),
				cvScalar(0, 0, 255));
	}
}

void Visualizer::getSliceCoor(CvPoint p, int *s, CvPoint *c) {
	*s = -1;
	*c = cvPoint(0, 0);
	if (outputImage) {
		int x, y, t, l;
		x = p.x / gw;
		y = p.y / gh;
		l = (p.x % gw) * width / gw;
		t = (p.y % gh) * height / gh;
		*s = startSlice + y * gx + x;
		if (*s <= SLICE_END && *s < startSlice + numSlices) {
			*c = cvPoint(l, t);
		} else {
			*s = -1;
			*c = cvPoint(0, 0);
		}
	}
}
