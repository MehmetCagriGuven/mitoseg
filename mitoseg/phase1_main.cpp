/*
 * phase1_main.cpp
 *
 *  Created on: Jan 13, 2017
 *      Author: fst
 */

#include "phase1_main.h"

#include <stdio.h>
#include <pthread.h>
#include <semaphore.h>

sem_t sem1;
pthread_mutex_t roilock = PTHREAD_MUTEX_INITIALIZER;

void mainFunctionPhase1(int s) {
	IplImage *originalImage = 0, *roiImage = 0;
	IplImage *grayImage = 0, *alImage = 0, *resmpImage = 0, *visualResmpImage =
			0, *visualResmpImage_color = 0;
	IplImage *bilateral = 0, *smoothed = 0, *visualBilateral = 0,
			*visualSmoothed = 0;
	ridge *r = 0;
	IplImage *ridges = 0, *visualRidges = 0, *ridgeMask = 0;
	IplImage *ridges_dir = 0, *visualRidges_dir = 0;
	localEnergy *e_lo = 0, *e_hi = 0;
	IplImage *visualLE_lo = 0, *visualLE_hi = 0;
	int lew_lo, leh_lo, lew_hi, leh_hi;
	curve *clist_lo = 0, *clist_hi = 0;
	curve *clist_lo_flt = 0, *clist_hi_flt = 0, *clist_hix = 0;
	int nc_lo, nc_hi, nc_lo_flt, nc_hi_flt, nc_hix;
	IplImage *visualCurves_lo = 0, *visualCurves_hi = 0;
	IplImage *visualCurves_hilo = 0;

	// Load Image Slice
	loadSlice(s, &originalImage);

	// Get image roi
	if (s != SLICE_START)
		pthread_mutex_lock(&roilock);
	getImageROI(originalImage, &roiImage); // roi image is created inside function as type of original image
	pthread_mutex_unlock(&roilock);

	// Convert to gray image
	grayImage = cvCreateImage(cvSize(roiImage->width, roiImage->height), 8, 1);
	cvCvtColor(roiImage, grayImage, CV_BGR2GRAY);

	// Auto levels
	alImage = cvCreateImage(cvSize(roiImage->width, roiImage->height), 8, 1);
	autoLevels(grayImage, alImage);

	// Resample image
	getResampledImage(alImage, &resmpImage); // resampled image is created inside function as IPL_DEPTH_32F

	// Create Core Images
	bilateral = cvCreateImage(cvSize(resmpImage->width, resmpImage->height),
	IPL_DEPTH_32F, 1);
	smoothed = cvCreateImage(cvSize(resmpImage->width, resmpImage->height),
	IPL_DEPTH_32F, 1);
	r = (ridge *) malloc(
			sizeof(ridge) * resmpImage->width * resmpImage->height);
	ridges = cvCreateImage(cvSize(resmpImage->width, resmpImage->height),
	IPL_DEPTH_32F, 1);
	ridges_dir = cvCreateImage(cvSize(resmpImage->width, resmpImage->height),
	IPL_DEPTH_32F, 1);
	ridgeMask = cvCreateImage(cvSize(resmpImage->width, resmpImage->height), 8,
			1);

	// Create Visual Images
	visualResmpImage = cvCreateImage(
			cvSize(resmpImage->width, resmpImage->height), 8, 1);
	visualResmpImage_color = cvCreateImage(
			cvSize(resmpImage->width, resmpImage->height), 8, 3);
	visualBilateral = cvCreateImage(
			cvSize(resmpImage->width, resmpImage->height), 8, 1);
	visualSmoothed = cvCreateImage(
			cvSize(resmpImage->width, resmpImage->height), 8, 1);
	visualRidges = cvCreateImage(cvSize(resmpImage->width, resmpImage->height),
			8, 1);
	visualRidges_dir = cvCreateImage(
			cvSize(resmpImage->width, resmpImage->height), 8, 3);
	//visualLE_lo = created inside function...
	//visualLE_hi = created inside function...
	visualCurves_lo = cvCreateImage(
			cvSize(resmpImage->width, resmpImage->height), 8, 3);
	visualCurves_hi = cvCreateImage(
			cvSize(resmpImage->width, resmpImage->height), 8, 3);
	visualCurves_hilo = cvCreateImage(
			cvSize(resmpImage->width, resmpImage->height), 8, 3);

	// Pre-processing
	// Smooth
	fastBilateralFilter(resmpImage, bilateral, (int) BSMOOTH_SPAT, BSMOOTH_GRAY,
			0);
	cvSmooth(bilateral, smoothed, CV_GAUSSIAN, 0, 0, SMOOTH_STDEV, 0);
	//ByPass:cvSmooth(resmpImage, smoothed, CV_GAUSSIAN, 0, 0, SMOOTH_STDEV, 0 );	//by-pass bilateral filtering
	/////////////////

	// Get Ridges and normalize (autolevels algorithm)
	getRidges(smoothed, r, ridges, ridges_dir, NULL, NULL, ridgeMask);
	normalizeRidges(r, ridges, ridgeMask);

	// Get local energy
	getLocalEnergyMap(r, ridges->width, ridges->height, (int) LE_SSIZE_LO,
			(int) LE_BSIZE_LO, &e_lo, &lew_lo, &leh_lo);
	getLocalEnergyMap(r, ridges->width, ridges->height, (int) LE_SSIZE_HI,
			(int) LE_BSIZE_HI, &e_hi, &lew_hi, &leh_hi);

	// Fit curves
	initCurveFitting();
	nc_lo = fitCurves(e_lo, lew_lo, leh_lo, &clist_lo, FCS_MAXCURVE_LO,
	FCS_XYSTEP_LO, FCS_XYRANGE_LO, FCS_INIT_THRESH_LO, FCS_MINLEN_LO);
	nc_hi = fitCurves(e_hi, lew_hi, leh_hi, &clist_hi, FCS_MAXCURVE_HI,
	FCS_XYSTEP_HI, FCS_XYRANGE_HI, FCS_INIT_THRESH_HI, FCS_MINLEN_HI);

	// Filter curves
	filterCurves(clist_lo, nc_lo, &clist_lo_flt, &nc_lo_flt,
			(float) FLC_THRESH_LO, (float) FLC_AVGTHRESH_LO);
	filterCurves(clist_hi, nc_hi, &clist_hi_flt, &nc_hi_flt,
			(float) FLC_THRESH_HI, (float) FLC_AVGTHRESH_HI);
	getHixCurves(clist_lo_flt, nc_lo_flt, lew_lo, leh_lo, clist_hi_flt,
			nc_hi_flt, &clist_hix, &nc_hix);

	// Visualize
	cvNormalize(resmpImage, visualResmpImage, 255, 0, CV_MINMAX);
	cvMerge(visualResmpImage, visualResmpImage, visualResmpImage, NULL,
			visualResmpImage_color);
	cvNormalize(bilateral, visualBilateral, 255, 0, CV_MINMAX);
	cvNormalize(smoothed, visualSmoothed, 255, 0, CV_MINMAX);
	cvNormalize(ridges, visualRidges, 255, 0, CV_MINMAX);
	visualizeRidgesDir(4, r, visualRidges_dir);
	visualizeLEMajority(e_lo, &visualLE_lo, lew_lo, leh_lo, VISUAL_LE_SSIZE, 1);
	visualizeLEMajority(e_hi, &visualLE_hi, lew_hi, leh_hi, VISUAL_LE_SSIZE, 1);
	visualizeCurves(visualResmpImage, clist_lo, nc_lo, (int) LE_SSIZE_LO,
			visualCurves_lo, (float) FLC_THRESH_LO, (float) FLC_AVGTHRESH_LO);
	visualizeCurves(visualResmpImage, clist_hi, nc_hi, (int) LE_SSIZE_HI,
			visualCurves_hi, (float) FLC_THRESH_HI, (float) FLC_AVGTHRESH_HI);
	visualizeHiLoCurves(clist_lo_flt, nc_lo_flt, lew_lo, leh_lo, clist_hi_flt,
			nc_hi_flt, CURV_HILO_COV_TH, visualResmpImage, visualCurves_hilo);

	// Save data
	saveImageDat(s, "resmpImage_", resmpImage);
	saveRidgeStructure(s, "ridgeStruc_", r,
			resmpImage->width * resmpImage->height);
	saveLocalEnergyStructure(s, "e_lo_", e_lo, lew_lo, leh_lo);
	saveLocalEnergyStructure(s, "e_hi_", e_hi, lew_hi, leh_hi);
	saveCurveStructure(s, "clist_lo_", clist_lo, nc_lo);
	saveCurveStructure(s, "clist_hi_", clist_hi, nc_hi);
	saveCurveStructure(s, "clist_lo_flt_", clist_lo_flt, nc_lo_flt);
	saveCurveStructure(s, "clist_hi_flt_", clist_hi_flt, nc_hi_flt);
	saveCurveStructure(s, "clist_hix_", clist_hix, nc_hix);

	// Save Visual Images
	saveImage(s, "00_image_", "", visualResmpImage);
	saveImage(s, "01_bilateral_", "", visualBilateral);
	saveImage(s, "02_smoothed_", "", visualSmoothed);
	saveImage(s, "03_ridges_", "", visualRidges);
	saveImage(s, "04_ridges_dir_", "", visualRidges_dir);
	saveImage(s, "05_le_lo_", "", visualLE_lo);
	saveImage(s, "06_le_hi_", "", visualLE_hi);
	saveImage(s, "07_curves_lo_", "", visualCurves_lo);
	saveImage(s, "08_curves_hi_", "", visualCurves_hi);
	saveImage(s, "09_curves_hilo_", "", visualCurves_hilo);

	// Free mem.
	if (originalImage)
		cvReleaseImage(&originalImage);
	if (roiImage)
		cvReleaseImage(&roiImage);
	if (grayImage)
		cvReleaseImage(&grayImage);
	if (alImage)
		cvReleaseImage(&alImage);
	if (resmpImage)
		cvReleaseImage(&resmpImage);
	if (smoothed)
		cvReleaseImage(&smoothed);
	if (ridges)
		cvReleaseImage(&ridges);
	if (ridges_dir)
		cvReleaseImage(&ridges_dir);
	if (ridgeMask)
		cvReleaseImage(&ridgeMask);
	if (visualResmpImage)
		cvReleaseImage(&visualResmpImage);
	if (visualResmpImage_color)
		cvReleaseImage(&visualResmpImage_color);
	if (visualSmoothed)
		cvReleaseImage(&visualSmoothed);
	if (visualRidges)
		cvReleaseImage(&visualRidges);
	if (visualRidges_dir)
		cvReleaseImage(&visualRidges_dir);
	if (visualLE_lo)
		cvReleaseImage(&visualLE_lo);
	if (visualLE_hi)
		cvReleaseImage(&visualLE_hi);
	if (visualCurves_lo)
		cvReleaseImage(&visualCurves_lo);
	if (visualCurves_hi)
		cvReleaseImage(&visualCurves_hi);
	if (visualCurves_hilo)
		cvReleaseImage(&visualCurves_hilo);
	//
	if (r)
		free(r);
	if (e_lo)
		free(e_lo);
	if (e_hi)
		free(e_hi);
	if (clist_lo)
		free(clist_lo);
	if (clist_hi)
		free(clist_hi);
	if (clist_lo_flt)
		free(clist_lo_flt);
	if (clist_hi_flt)
		free(clist_hi_flt);
	if (clist_hix)
		free(clist_hix);
}

void *threadPhase1(void *t_data) {
	int s = *(int *) t_data;

	printf("Processing slice #%d...\n", s);
	mainFunctionPhase1(s);

	sem_post(&sem1);
	pthread_exit(NULL);
}

void startPhase1() {
	int s, t;
	int rc;
	int numThreads = SLICE_END - SLICE_START + 1;
	pthread_t *threads = (pthread_t *) malloc(sizeof(pthread_t) * numThreads);
	pthread_attr_t attr;
	int *t_data = (int *) malloc(sizeof(int) * numThreads);
	sem_init(&sem1, 0, numCores);
	pthread_mutex_trylock(&roilock);
	void *status;

	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

	printf(">>>> PHASE #1\n");
	for (s = SLICE_START, t = 0; s <= SLICE_END; s++, t++) {
		sem_wait(&sem1);
		t_data[t] = s;
		rc = pthread_create(&threads[t], &attr, threadPhase1,
				(void *) &t_data[t]);
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

	free(threads);
	free(t_data);
}
