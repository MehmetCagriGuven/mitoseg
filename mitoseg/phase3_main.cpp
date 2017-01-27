/*
 * phase3_main.cpp
 *
 *  Created on: Jan 15, 2017
 *      Author: fst
 */

#include "phase3_main.h"

void revalidateSnakes(snake25d *snakeArray, int n) {
	// Visualizer
	Visualizer v;
	v.setStartSlice(SLICE_START);
	v.setNumSlices(sn25d_t);
	v.setFileNameTag(DESTPATH, "00_image_", "");
	v.update();

	// w, h
	int w_lo, h_lo;
	int w_hi, h_hi;
	getMapSize(v.getFrameSize(), (int) LE_SSIZE_LO, &w_lo, &h_lo);
	getMapSize(v.getFrameSize(), (int) LE_SSIZE_HI, &w_hi, &h_hi);

	// Load & process curve data
	dataPacket dp;
	dp.cs_lo = loadCurveStack(SLICE_START, SLICE_END, "clist_lo_flt_");
	dp.cs_hi = loadCurveStack(SLICE_START, SLICE_END, "clist_hix_");
	dp.ces_lo = createCurveEnergyStack(dp.cs_lo, w_lo, h_lo);
	dp.ces_hi = createCurveEnergyStack(dp.cs_hi, w_hi, h_hi);
	dp.grad_ces_lo = createCurveEnergyGradientStack(dp.ces_lo);
	dp.cps_lo = createCurvePointsStack(dp.cs_lo);
	dp.cps_hi = createCurvePointsStack(dp.cs_hi, w_hi, h_hi);

	for (int i = 0; i < n; i++) {
		validateSnake25d(&dp, &snakeArray[i]);
		printf("\rRevalidating... %d / %d ", i + 1, n);
	}
	printf("\n");

	// Free mem
	destroyCurvePointsStack(dp.cps_hi);
	destroyCurveEnergyGradientStack(dp.grad_ces_lo);
	destroyCurveEnergyStack(dp.ces_lo);
	destroyCurveEnergyStack(dp.ces_hi);
	destroyCurveStack(dp.cs_lo);
	destroyCurveStack(dp.cs_hi);
}

void startPhase3() {
	int start_z = 0, end_z = 0, s;
	Visualizer v;

	printf(">>>> PHASE #3\n");

	// Visualize initial points
	CvPoint *pts = NULL;
	int n_pts;
	int t_n_pts = 0;
	for (s = SLICE_START; s + sn25d_t - 1 <= SLICE_END; s += sn25d_t)
	{
		start_z = s;
		end_z = s + sn25d_t - 1;
		// Load initial pts
		n_pts = loadPoints(start_z, end_z, &pts);
		t_n_pts += (n_pts * sn25d_t);
		// Visualizer
		v.setStartSlice(start_z);
		v.setNumSlices(1);
		v.setFileNameTag(DESTPATH, "00_image_", "");
		v.setMarkerBuf(pts, n_pts, start_z, LE_SSIZE_LO);
		v.update();
		saveImage(start_z, "10_initPts_", "", v.outputImage);
	}
	v.setMarkerBuf(NULL);
	if (pts)
		free(pts);

	// Load snakes
	int z1 = SLICE_START;
	int z2 = end_z;
	snake25d *snakeArray = NULL;
	int n_snk = 0;
	poly25d *polyArray = NULL;
	int n_poly = 0;
	poly25d *polyArrayMerged = NULL;
	int n_poly_merged = 0;
	n_snk = loadSnakeArray(z1, z2, &snakeArray);
	printf("%d snakes loaded...\n", n_snk);
	//revalidateSnakes(snakeArray, n_snk);

	// Get valid polygons
	n_poly = convertValidSnakesToPolyArray(snakeArray, n_snk,
			POLY_VALIDITY - EPS, &polyArray);
	printf("%d valid polygons extracted...\n", n_poly);

	// Get merged polygons
	polyArrayMerged = (poly25d *) calloc(n_poly, sizeof(poly25d));
	n_poly_merged = mergeArrayOfPoly25d(polyArray, n_poly, POLY_MERGE,
			polyArrayMerged);
	printf("%d merged polygons extracted...\n", n_poly_merged);

	// Visualize polygons
	v.setFileNameTag(DESTPATH, "00_image_", "");
	v.setPolyBuf(polyArray, n_poly, LE_SSIZE_LO);
	for (s = z1; s <= z2; s++) {
		v.setStartSlice(s);
		v.setNumSlices(1);
		v.update();
		saveImage(s, "11_mitos_", "", v.outputImage);
	}
	v.setPolyBuf(polyArrayMerged, n_poly_merged, LE_SSIZE_LO);
	for (s = z1; s <= z2; s++) {
		v.setStartSlice(s);
		v.setNumSlices(1);
		v.update();
		saveImage(s, "12_mitos_merged_", "", v.outputImage);
	}
	savePolyArrayAsPLY(polyArrayMerged, n_poly_merged);

	// Free mem.
	if (snakeArray) {
		free(snakeArray);
		snakeArray = NULL;
	}
	if (polyArray && n_poly) {
		deinitPoly25dArray(polyArray, n_poly);
		free(polyArray);
		polyArray = NULL;
		n_poly = 0;
	}
	if (polyArrayMerged && n_poly_merged) {
		deinitPoly25dArray(polyArrayMerged, n_poly_merged);
		free(polyArrayMerged);
		polyArrayMerged = NULL;
		n_poly_merged = 0;
	}
}
