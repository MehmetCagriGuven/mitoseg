/*
 * phase2_main.cpp
 *
 *  Created on: Jan 15, 2017
 *      Author: fst
 */

#include "phase2_main.h"

void startPhase2() {
	printf(">>>> PHASE #2\n");

	// Initial Points
	CvPoint *pts = NULL;
	int n_pts;

	// Snakes
	int n_snk = 0;
	snake25dList *snakeList = NULL;
	snake25d *snakeArray = NULL;

	// w, h
	int w_lo, h_lo;
	int w_hi, h_hi;
	IplImage *rsmpImg = NULL;
	loadImageDat(SLICE_START, "resmpImage_", &rsmpImg);
	getMapSize(cvSize(rsmpImg->width, rsmpImg->height), (int) LE_SSIZE_LO,
			&w_lo, &h_lo);
	getMapSize(cvSize(rsmpImg->width, rsmpImg->height), (int) LE_SSIZE_HI,
			&w_hi, &h_hi);
	cvReleaseImage(&rsmpImg);

	// Load & process curve data
	dataPacket dp;
	dp.cs_lo = loadCurveStack(SLICE_START, SLICE_END, "clist_lo_flt_");
	dp.cs_hi = loadCurveStack(SLICE_START, SLICE_END, "clist_hix_");
	dp.ces_lo = createCurveEnergyStack(dp.cs_lo, w_lo, h_lo);
	dp.ces_hi = createCurveEnergyStack(dp.cs_hi, w_hi, h_hi);
	dp.grad_ces_lo = createCurveEnergyGradientStack(dp.ces_lo);
	dp.cps_lo = createCurvePointsStack(dp.cs_lo);
	dp.cps_hi = createCurvePointsStack(dp.cs_hi, w_hi, h_hi);

	int s, start_z, end_z = 0;

	for (s = SLICE_START; s + sn25d_t - 1 <= SLICE_END; s += sn25d_t)
	{
		start_z = s;
		end_z = s + sn25d_t - 1;

		n_pts = findInitialPointsWithinZRange(&dp, &pts, start_z, end_z);

		// Save data
		printf("Saving %d initial points...\n", n_pts);
		savePoints(pts, n_pts, start_z, end_z);

		n_snk += retrieveAllSnakes25d(&dp, pts, n_pts, start_z, end_z,
				&snakeList);
	}

	// Save data
	convertSnakeList2Array(snakeList, &snakeArray);
	printf("Saving %d snakes...\n", n_snk);
	saveSnakeArray(snakeArray, n_snk, SLICE_START, end_z);

	// Free mem
	destroyCurvePointsStack(dp.cps_hi);
	destroyCurveEnergyGradientStack(dp.grad_ces_lo);
	destroyCurveEnergyStack(dp.ces_lo);
	destroyCurveEnergyStack(dp.ces_hi);
	destroyCurveStack(dp.cs_lo);
	destroyCurveStack(dp.cs_hi);
	destroySnakeList(&snakeList);
	if (pts)
		free(pts);
	if (snakeArray)
		free(snakeArray);
}
