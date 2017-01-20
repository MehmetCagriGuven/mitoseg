/*
 * preprocess.cpp
 *
 *  Created on: Jan 13, 2017
 *      Author: fst
 */

#include "preprocess.h"

void autoLevels(IplImage *input, IplImage *output) {
	int i, a, b, x, y;
	uchar *row, *col, *orow, *ocol;
	double prob[256], p, c, m;
	int bins[] = { 258 };
	float range[] = { -1, 256 };
	float* ranges[] = { range };
	CvHistogram* hist;
	hist = cvCreateHist(1, bins, CV_HIST_ARRAY, ranges, 1);
	IplImage* im[] = { input };

	cvCalcHist(im, hist);
	for (i = 0; i < 256; i++)
		prob[i] = (double) cvQueryHistValue_1D(hist, i + 1)
				/ (input->width * input->height);
	a = 0;
	p = 0;
	while (p < AUTOLEVELSCUT) {
		p += prob[a];
		a++;
	}
	b = 255;
	p = 0;
	while (p < AUTOLEVELSCUT) {
		p += prob[b];
		b--;
	}
	m = 255.0 / (b - a);
	c = -m * a;
	for (y = 0, row = (uchar *) input->imageData, orow =
			(uchar *) output->imageData; y < input->height;
			y++, row += input->widthStep, orow += output->widthStep) {
		for (x = 0, col = row, ocol = orow; x < input->width;
				x++, col++, ocol++) {
			i = cvRound(*col * m + c);
			if (i < 0)
				i = 0;
			else if (i > 255)
				i = 255;
			*ocol = i;
		}
	}
	cvReleaseHist(&hist);
}

void getResampledImage(IplImage *input, IplImage **output) {
	double sf;
	IplImage *temp;
	if (*output)
		cvReleaseImage(output);
	if (TFACTOR == RESOLUTION) {
		*output = cvCreateImage(cvSize(input->width, input->height),
		IPL_DEPTH_32F, 1);
		cvConvert(input, *output);
	} else {
		sf = RESOLUTION / TFACTOR;
		temp = cvCreateImage(cvSize(input->width, input->height), IPL_DEPTH_32F,
				1);
		*output = cvCreateImage(
				cvSize(cvRound(input->width * sf), cvRound(input->height * sf)),
				IPL_DEPTH_32F, 1);
		cvConvert(input, temp);
		cvResize(temp, *output, CV_INTER_LINEAR);
		cvReleaseImage(&temp);
	}
}

void getRidges(IplImage *input, ridge *r, IplImage *output,
		IplImage *output_dir, IplImage *output_dir4, IplImage *output_dir8,
		IplImage *mask) {
	IplImage *image = 0;
	IplImage *gxx = 0, *gxy = 0, *gyx = 0, *gyy = 0;
	CvMat *V, *E;

	unsigned int row, ptr, mrow, mptr;
	int i, j, ind;
	float m[4], p, a;
	CvMat M;

	image = cvCreateImage(cvSize(input->width, input->height), IPL_DEPTH_32F,
			1);
	gxx = cvCreateImage(cvSize(input->width, input->height), IPL_DEPTH_32F, 1);
	gxy = cvCreateImage(cvSize(input->width, input->height), IPL_DEPTH_32F, 1);
	gyx = cvCreateImage(cvSize(input->width, input->height), IPL_DEPTH_32F, 1);
	gyy = cvCreateImage(cvSize(input->width, input->height), IPL_DEPTH_32F, 1);
	V = cvCreateMat(2, 2, CV_32F);
	E = cvCreateMat(2, 1, CV_32F);

	cvCopy(input, image);

	// Ridge detection using Hessian
	cvSobel(image, gxx, 2, 0, 3);
	cvSobel(image, gxy, 1, 1, 3);
	cvSobel(image, gyx, 1, 1, 3);
	cvSobel(image, gyy, 0, 2, 3);

	for (i = 0, row = 0, mrow = 0, ind = 0; i < image->height;
			i++, row += image->widthStep, mrow += mask->widthStep) {
		for (j = 0, ptr = row, mptr = mrow; j < image->width; j++, ptr +=
				sizeof(float), mptr++, ind++) {
			m[0] = *(float *) (gxx->imageData + ptr);
			m[1] = *(float *) (gxy->imageData + ptr);
			m[2] = *(float *) (gyx->imageData + ptr);
			m[3] = *(float *) (gyy->imageData + ptr);

			M = cvMat(2, 2, CV_32F, m);
			cvEigenVV(&M, V, E);

			if (E->data.fl[0] < 0) {
				p = 0;
				*(char *) (mask->imageData + mptr) = (char) 0;
				r[ind].mask = 0;
			} else {
				if (E->data.fl[1] < 0)
					p = E->data.fl[0];
				else
					p = E->data.fl[0] - E->data.fl[1];
				*(char *) (mask->imageData + mptr) = (char) 255;
				r[ind].mask = (char) 255;
			}

			r[ind].mag = p;

			a = (float) (atan2(-V->data.fl[1], V->data.fl[0]) * 180 / PI);
			if (a < 0)
				a += 360;
			r[ind].angle = a;
			r[ind].vec[0] = V->data.fl[0];
			r[ind].vec[1] = V->data.fl[1];
			r[ind].val[0] = E->data.fl[0];
			r[ind].val[1] = E->data.fl[1];
			r[ind].dir = (a > 180) ? a - 180 : a;

			// dir4
			if ((a > 22.5 && a <= 67.5) || (a > 202.5 && a <= 247.5))
				r[ind].dir4 = 1;
			else if ((a > 67.5 && a <= 112.5) || (a > 247.5 && a <= 292.5))
				r[ind].dir4 = 2;
			else if ((a > 112.5 && a <= 157.5) || (a > 292.5 && a <= 337.5))
				r[ind].dir4 = 3;
			else
				r[ind].dir4 = 0;

			// dir8
			if ((a > 11.25 && a <= 33.75) || (a > 191.25 && a <= 213.75))
				r[ind].dir8 = 1;
			else if ((a > 33.75 && a <= 56.25) || (a > 213.75 && a <= 236.25))
				r[ind].dir8 = 2;
			else if ((a > 56.25 && a <= 78.75) || (a > 236.25 && a <= 258.75))
				r[ind].dir8 = 3;
			else if ((a > 78.75 && a <= 101.25) || (a > 258.75 && a <= 281.25))
				r[ind].dir8 = 4;
			else if ((a > 101.25 && a <= 123.75) || (a > 281.25 && a <= 303.75))
				r[ind].dir8 = 5;
			else if ((a > 123.75 && a <= 146.25) || (a > 303.75 && a <= 326.25))
				r[ind].dir8 = 6;
			else if ((a > 146.25 && a <= 168.75) || (a > 326.25 && a <= 348.75))
				r[ind].dir8 = 7;
			else
				r[ind].dir8 = 0;

			if (output_dir8) {
				cvSetReal2D(output_dir8, i, j, r[ind].dir8);
			}

			if (output_dir4) {
				cvSetReal2D(output_dir4, i, j, r[ind].dir4);
			}

			if (output_dir) {
				cvSetReal2D(output_dir, i, j, r[ind].dir);
			}

			cvSetReal2D(output, i, j, r[ind].mag);
		}
	}

	cvReleaseImage(&image);
	cvReleaseImage(&gxx);
	cvReleaseImage(&gxy);
	cvReleaseImage(&gyx);
	cvReleaseImage(&gyy);
	cvReleaseMat(&E);
	cvReleaseMat(&V);
}

void visualizeRidgesDir(int opt, ridge *r, IplImage *visual_dir) {
	int i, j, ind, mind = visual_dir->width * visual_dir->height, p;
	float ridges_max = r[0].mag, ridges_min = r[0].mag;

	for (ind = 0; ind < mind; ind++) {
		if (ridges_max < r[ind].mag)
			ridges_max = r[ind].mag;
		if (ridges_min > r[ind].mag)
			ridges_min = r[ind].mag;
	}

	for (i = 0, ind = 0; i < visual_dir->height; i++) {
		for (j = 0; j < visual_dir->width; j++, ind++) {
			p = cvRound(
					255.0f * (r[ind].mag - ridges_min)
							/ (ridges_max - ridges_min));

			if (opt == 8) {
				if (r[ind].dir8 == 7)
					cvSet2D(visual_dir, i, j, cvScalar(p, 0, 0));
				else if (r[ind].dir8 == 6)
					cvSet2D(visual_dir, i, j, cvScalar(0, p, 0));
				else if (r[ind].dir8 == 5)
					cvSet2D(visual_dir, i, j, cvScalar(0, 0, p));
				else if (r[ind].dir8 == 4)
					cvSet2D(visual_dir, i, j, cvScalar(p, p, 0));
				else if (r[ind].dir8 == 3)
					cvSet2D(visual_dir, i, j, cvScalar(p, 0, p));
				else if (r[ind].dir8 == 2)
					cvSet2D(visual_dir, i, j, cvScalar(0, p, p));
				else if (r[ind].dir8 == 1)
					cvSet2D(visual_dir, i, j, cvScalar(p, p, p));
				else
					cvSet2D(visual_dir, i, j, cvScalar(p, p, p / 2));
			} else if (opt == 4) {
				if (r[ind].dir4 == 3)
					cvSet2D(visual_dir, i, j, cvScalar(p, 0, 0));
				else if (r[ind].dir4 == 2)
					cvSet2D(visual_dir, i, j, cvScalar(0, p, 0));
				else if (r[ind].dir4 == 1)
					cvSet2D(visual_dir, i, j, cvScalar(0, 0, p));
				else
					cvSet2D(visual_dir, i, j, cvScalar(p, p, p));
			}
		}
	}
}

void normalizeRidges(ridge *r, IplImage *ridges, IplImage *mask) {
	const double autoLevelCut = 2.81;	// 0.5% cut for norm.dist.

	CvScalar mean, std;
	float p;
	int x, y, i;
	uchar *row, *col;
	cvAvgSdv(ridges, &mean, &std, mask);
	double a, b, m, c;

	a = mean.val[0] - std.val[0] * autoLevelCut;
	if (a < 0)
		a = 0;
	b = mean.val[0] + std.val[0] * autoLevelCut;
	m = 1. / (b - a);
	c = -m * a;

	for (y = 0, i = 0, row = (uchar *) ridges->imageData; y < ridges->height;
			y++, row += ridges->widthStep) {
		for (x = 0, col = row; x < ridges->width;
				x++, i++, col += sizeof(float)) {
			if (r[i].mask) {
				p = (float) (m * r[i].mag + c);
				if (p > 1.0f)
					p = 1.0f;
				else if (p < 0)
					p = 0;

				r[i].mag = p;
				*((float *) col) = p;
			}

		}

	}

}

void fastBilateralFilter(IplImage *input, IplImage *output, int spatialWindow,
		float graySigma, int order) {
	const int bins = 64;

	spatialWindow |= 1;
	graySigma *= bins;

	int i, j, o, k;
	int imW, imH, imin, imax, jmin, jmax, histH, histWB, histHWB;
	int hk = spatialWindow / 2;

	double minval, maxval, valscale;
	double sigma_2, c;
	float grayk[bins];
	float conv, norm, term;
	int line_hist[bins];
	int diff;
	int *im, *im_imin_jmin, *im_imin_jmax, *im_imax_jmin, *im_imax_jmax,
			*im_i_jmax, *im_i_j, *im_ii_j, *im_i_jj, *im_ii_jj, *im_imin_j,
			*im_imax_j, *im_i_jmin, *im_i, *im_ci, *im_ci_cj;
	int *ihist, *ihist_h, *ihist_h_1, *ihist_wh, *ihist_h_j, *ihist_h_1_j,
			*ihist_wh_j, *ihist_h_wj, *ihist_wh_wj, *ihist_end;
	int hist;
	uchar *row, *col;

	// Get value range & statistics
	cvMinMaxLoc(input, &minval, &maxval);
	valscale = (bins - 1) / (maxval - minval);

	// Init Gaussian gray kernel
	c = 1.0 / (sqrt(2.0 * PI) * graySigma);
	sigma_2 = 2.0 * graySigma * graySigma;
	for (i = 0; i < bins; i++) {
		grayk[i] = (float) (c * exp(-(double) (i * i) / sigma_2));
	}

	// Init integral image & histogram
	imW = input->width + spatialWindow + 1;
	imH = input->height + spatialWindow + 1;
	histH = spatialWindow + 1;
	histWB = imW * bins;
	histHWB = histH * histWB;
	jmin = imin = hk + 1;
	jmax = imW - hk - 2;
	imax = imH - hk - 2;
	im = (int *) malloc(sizeof(int) * imH * imW);
	memset((void *) im, 0, sizeof(int) * imH * imW);
	ihist = (int *) malloc(sizeof(int) * histH * histWB);
	memset((void *) ihist, 0, sizeof(int) * histH * histWB);
	im_imin_jmin = im + imin * imW + jmin;
	im_imin_jmax = im + imin * imW + jmax;
	im_imax_jmin = im + imax * imW + jmin;
	im_imax_jmax = im + imax * imW + jmax;
	ihist_end = ihist + histHWB;

	// Quantize image levels to integers
	for (row = (uchar *) input->imageData, im_i_jmin = im_imin_jmin;
			im_i_jmin <= im_imax_jmin;
			im_i_jmin += imW, row += input->widthStep) {
		im_i_jmax = im_i_jmin + input->width - 1;
		for (col = row, im_i_j = im_i_jmin; im_i_j <= im_i_jmax;
				im_i_j++, col += sizeof(float)) {
			*im_i_j = (int) (((double) (*(float *) col) - minval) * valscale);
		}
	}

	// Order loop (box, triangle, approx. Gaussian) filters
	for (o = 0; o <= order; o++) {
		// Prepare image borders
		for (j = jmin, im_imin_j = im_imin_jmin, im_imax_j = im_imax_jmin;
				j <= jmax; j++, im_imin_j++, im_imax_j++) {
			for (i = 1, im_i_j = im_imin_j - imW, im_ii_j = im_imax_j + imW;
					i < imin; i++, im_i_j -= imW, im_ii_j += imW) {
				*im_i_j = *im_imin_j;
				*im_ii_j = *im_imax_j;
			}
		}
		for (i = imin, im_i_jmin = im_imin_jmin, im_i_jmax = im_imin_jmax;
				i <= imax; i++, im_i_jmin += imW, im_i_jmax += imW) {
			for (j = 1, im_i_j = im_i_jmin - 1, im_i_jj = im_i_jmax + 1;
					j < jmin; j++, im_i_j--, im_i_jj++) {
				*im_i_j = *im_i_jmin;
				*im_i_jj = *im_i_jmax;
			}
		}
		for (i = 1, im_i_j = im_imin_jmin - imW, im_ii_j = im_imax_jmin + imW, im_i_jj =
				im_imin_jmax - imW, im_ii_jj = im_imax_jmax + imW; i < imin;
				i++, im_i_j -= imW, im_ii_j += imW, im_i_jj -= imW, im_ii_jj +=
						imW) {
			for (j = 1; j < jmin; j++) {
				*(im_i_j - j) = *im_imin_jmin;
				*(im_ii_j - j) = *im_imax_jmin;
				*(im_i_jj + j) = *im_imin_jmax;
				*(im_ii_jj + j) = *im_imax_jmax;
			}
		}

		// Init pointers
		im_i = im + imW;
		im_ci = im + (1 - hk) * imW;
		ihist_h = ihist + histWB;
		ihist_h_1 = ihist;
		ihist_wh = ihist + 2 * histWB;
		row = (uchar *) output->imageData
				+ (1 - spatialWindow) * input->widthStep;

		// Filter loop
		for (i = 1; i < imH - 1; i++) {
			// Clear temp histogram item
			memset((void *) line_hist, 0, sizeof(int) * bins);

			for (j = 1, im_i_j = im_i + 1, im_ci_cj = im_ci + 1 - hk, ihist_h_j =
					ihist_h + bins, ihist_h_1_j = ihist_h_1 + bins, ihist_wh_j =
					ihist_wh + bins, ihist_h_wj = ihist_h
					+ (1 - spatialWindow) * bins, ihist_wh_wj = ihist_wh
					+ (1 - spatialWindow) * bins, col = row
					+ (1 - spatialWindow) * sizeof(float); j < imW - 1;
					j++, im_i_j++, im_ci_cj++, ihist_h_j += bins, ihist_h_1_j +=
							bins, ihist_wh_j += bins, ihist_h_wj += bins, ihist_wh_wj +=
							bins, col += sizeof(float)) {
				// Calculate integral histogram
				line_hist[*im_i_j]++;
				for (k = 0; k < bins; k++)
					ihist_h_j[k] = line_hist[k] + ihist_h_1_j[k];

				// Calculate filtered image & hist
				if (i >= spatialWindow && j >= spatialWindow) {
					norm = 0;
					conv = 0;
					for (k = 0; k < bins; k++) {
						hist = ihist_h_j[k] - ihist_wh_j[k] - ihist_h_wj[k]
								+ ihist_wh_wj[k];
						diff = (*im_ci_cj > k) ? *im_ci_cj - k : k - *im_ci_cj;
						term = hist * grayk[diff];
						conv += term * k;
						norm += term;
					}
					*(float *) col = conv / norm;
				}
			}

			// Update pointers
			im_i += imW;
			im_ci += imW;
			ihist_h += histWB;
			if (ihist_h >= ihist_end)
				ihist_h -= histHWB;
			ihist_h_1 += histWB;
			if (ihist_h_1 >= ihist_end)
				ihist_h_1 -= histHWB;
			ihist_wh += histWB;
			if (ihist_wh >= ihist_end)
				ihist_wh -= histHWB;
			row += output->widthStep;
		}
	}

	// Free mem.
	free(im);
	free(ihist);
}

void saveRidgeStructure(int s, char const *tag, ridge *r, int numRecords) {
	FILE *f;
	char outputfn[256];
	char outputsrcfn[256];
	sprintf(outputsrcfn, "%s%s%s.dat", DESTPATH, tag, FNAME);
	sprintf(outputfn, outputsrcfn, s);
	if (r) {
		printf("Saving: %s\n", outputfn);
		f = fopen(outputfn, "wb");
		fwrite(&numRecords, sizeof(int), 1, f);
		fwrite(r, sizeof(ridge), numRecords, f);
		fclose(f);
	}
}

int loadRidgeStructure(int s, char const *tag, ridge **r) {
	int numRecords;
	FILE *f;
	size_t t;
	char outputfn[256];
	char outputsrcfn[256];
	sprintf(outputsrcfn, "%s%s%s.dat", DESTPATH, tag, FNAME);
	sprintf(outputfn, outputsrcfn, s);
	if (*r) {
		free(*r);
	}
	printf("Loading: %s\n", outputfn);
	f = fopen(outputfn, "rb");
	t = fread(&numRecords, sizeof(int), 1, f);
	*r = (ridge *) malloc(sizeof(ridge) * numRecords);
	t = fread(*r, sizeof(ridge), numRecords, f);
	if (!t)
		printf("Error loading data!\n");
	fclose(f);
	return numRecords;
}
