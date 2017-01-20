/*
 * curve.cpp
 *
 *  Created on: Jan 13, 2017
 *      Author: fst
 */

#include "curve.h"

double angle_x[LE_BINS];
double angle_y[LE_BINS];

void getLocalEnergyMap(ridge *r, const int width, const int height,
		const int ssize, const int bsize, localEnergy **e, int *w, int *h) {
	*w = (width / ssize) + ((width % ssize) ? 1 : 0);
	*h = (height / ssize) + ((height % ssize) ? 1 : 0);
	int maxi = (*w) * (*h);
	int maxj = width * height;
	int hsize = bsize / 2;
	int x, y, x1, y1, x2, y2, xx, yy, i, j, b, major;
	float m;

	// Init bin indices
	int *q = (int *) malloc(sizeof(int) * maxj);
	for (j = 0; j < maxj; j++) {
		q[j] = cvRound(LE_BINS * r[j].angle / 180.0f) % LE_BINS;
	}

	// Alloc mem.
	if (*e)
		free(*e);
	*e = (localEnergy *) malloc(sizeof(localEnergy) * maxi);
	memset(*e, 0, sizeof(localEnergy) * maxi);

	// Main loop
	for (y = 0, i = 0; y < height; y += ssize) {
		y1 = y - hsize;
		y2 = y + hsize;
		if (y1 < 0)
			y1 = 0;
		if (y2 > height)
			y2 = height;

		for (x = 0; x < width; x += ssize, i++) {
			x1 = x - hsize;
			x2 = x + hsize;
			if (x1 < 0)
				x1 = 0;
			if (x2 > width)
				x2 = width;

			// Calculate energy per direction
			for (yy = y1; yy < y2; yy++) {
				for (xx = x1; xx < x2; xx++) {
					j = yy * width + xx;
					if (r[j].mask) {
						(*e)[i].energy[q[j]] += r[j].mag;
					}
				}
			}

			// Find major direction
			major = 0;
			m = 0;
			for (b = 0; b < LE_BINS; b++) {
				if ((*e)[i].energy[b] > m) {
					m = (*e)[i].energy[b];
					major = b;
				}
			}
			(*e)[i].major = major;
		}
	}

	// Free mem.
	free(q);
}

void visualizeLEMajority(localEnergy *e, IplImage **visual, const int w,
		const int h, const int ssize, const int showDir) {
	float maxe = 0;
	double a;
	int i, x, y, b, k, xx, yy, mx, my, x1, x2, y1, y2, dls;

	if (*visual) {
		cvReleaseImage(visual);
	}
	*visual = cvCreateImage(cvSize(ssize * w, ssize * h), 8, 1);

	dls = ssize / 2 - 2;
	if (dls < 1)
		dls = 1;

	// Find global highest energy
	for (y = 0, i = 0; y < h; y++) {
		for (x = 0; x < w; x++, i++) {
			for (b = 0; b < LE_BINS; b++) {
				if (e[i].energy[b] > maxe)
					maxe = e[i].energy[b];
			}
		}
	}

	// Main loop
	for (y = 0, i = 0; y < h; y++) {
		yy = y * ssize;

		for (x = 0; x < w; x++, i++) {
			xx = x * ssize;
			k = cvRound(255 * e[i].energy[e[i].major] / maxe);
			cvRectangle(*visual, cvPoint(xx, yy),
					cvPoint(xx + ssize - 1, yy + ssize - 1), cvScalar(k),
					CV_FILLED);
			if (showDir) {
				a = (PI / 2) - (PI * e[i].major) / LE_BINS;
				mx = cvRound(cos(a) * dls);
				my = cvRound(sin(a) * dls);
				x1 = x2 = xx + ssize / 2;
				y1 = y2 = yy + ssize / 2;
				x1 += mx;
				y1 += my;
				x2 -= mx;
				y2 -= my;
				if (k > 128)
					k = 0;
				else
					k = 255;
				cvLine(*visual, cvPoint(x1, y1), cvPoint(x2, y2), cvScalar(k));
			}
		}
	}
}

void calcCurveParams(localEnergy *e, const int w, const int h, curve *c,
		const double S) {

	double t, k, dt, m, wtheta;
	double x, y, xp, yp, xn, yn, mx, my;
	int xi, yi, i, j;

	// Calc. model params
	int dx = c->x2 - c->x1;
	int dy = c->y2 - c->y1;
	double R = sqrt((double) (dx * dx + dy * dy));
	double b = (double) (4 * c->h) / R;
	double a = -b / R;
	double cost = dx / R;
	double sint = dy / R;

	c->R = R;
	c->a = a;
	c->b = b;
	c->sint = sint;
	c->cost = cost;
	c->score = 0;
	c->len = 0;
	if (R < EPS)
		return;

	// Main loop
	t = -S / sqrt(1.0 + b * b);
	k = a * t * t + b * t;
	xp = cost * t - sint * k + c->x1;
	yp = sint * t + cost * k + c->y1;
	x = (double) c->x1;
	y = (double) c->y1;
	dt = -t;
	t = 0;
	while (1) {
		t += dt;
		k = a * t * t + b * t;

		xn = cost * t - sint * k + c->x1;
		yn = sint * t + cost * k + c->y1;

		mx = xn - xp;
		my = yn - yp;

		xi = cvRound(x);
		yi = cvRound(y);
		if (xi >= 0 && yi >= 0 && xi < w && yi < h) {
			i = yi * w + xi;
			for (j = 0; j < LE_BINS; j++) {
				wtheta = angle_x[j] * my + angle_y[j] * mx;
				wtheta = 2 * (wtheta * wtheta) / (mx * mx + my * my) - 1;
				c->score += (float) wtheta * e[i].energy[j];
			}
		}

		c->len++;

		if (t > R)
			break;

		xp = x;
		yp = y;
		x = xn;
		y = yn;
		m = 2 * a * t + b;
		dt = S / sqrt(1.0 + m * m);
	}

}

void calcAuxCurveParams(curve *c) {
	c->avgscore = c->score / c->len;
}

void drawCurve(IplImage *im, curve *c, const int ssize, const int psize,
		const double S, CvScalar cl) {
	double t, k, dt, m, x, y;
	int hsize = ssize / 2;

	t = 0;
	while (t <= c->R && c->R >= EPS) {
		k = c->a * t * t + c->b * t;
		x = c->cost * t - c->sint * k + c->x1;
		y = c->sint * t + c->cost * k + c->y1;

		cvCircle(im,
				cvPoint(cvRound(x * ssize) + hsize, cvRound(y * ssize) + hsize),
				psize, cl, CV_FILLED);

		m = 2 * c->a * t + c->b;
		dt = S / sqrt(1.0 + m * m);
		t += dt;
	}
}

void drawCurve1(IplImage *im, curve *c, const int ssize, const int psize,
		const double S, CvScalar cl) {
	double t, k, dt, m, x, y;
	int hsize = 0;
	CvPoint p1, p2;
	int first = 1;

	t = 0;
	while (t <= c->R && c->R >= EPS) {
		k = c->a * t * t + c->b * t;
		x = c->cost * t - c->sint * k + c->x1;
		y = c->sint * t + c->cost * k + c->y1;

		if (first) {
			first = 0;
			p2 = cvPoint(cvRound(x * ssize) + hsize,
					cvRound(y * ssize) + hsize);
		} else {
			p1 = p2;
			p2 = cvPoint(cvRound(x * ssize) + hsize,
					cvRound(y * ssize) + hsize);
			cvLine(im, p1, p2, cl, 1);
		}

		m = 2 * c->a * t + c->b;
		dt = S / sqrt(1.0 + m * m);
		t += dt;
	}
}

void drawCurve2(IplImage *im, curve *c, const int ssize, const int psize,
		const double S, CvScalar cl) {
	double t, k, dt, m, x, y;
	int hsize = ssize / 2;

	t = -c->R;
	while (t <= 2 * c->R && c->R >= EPS) {
		k = c->a * t * t + c->b * t;
		x = c->cost * t - c->sint * k + c->x1;
		y = c->sint * t + c->cost * k + c->y1;

		cvCircle(im,
				cvPoint(cvRound(x * ssize) + hsize, cvRound(y * ssize) + hsize),
				psize, (0 <= t && t < c->R) ? cl : cvScalar(0, 0, 255),
				CV_FILLED);

		m = 2 * c->a * t + c->b;
		dt = S / sqrt(1.0 + m * m);
		t += dt;
	}
}

int fitCurve3(localEnergy *e, const int w, const int h, curve *c,
		const int maxIter, char *mask)	// Currently used version
		{
	int i, dx, dy;
	curve best, cur, max, n;

	const double maxcurv = FC_MAXCURV * FC_MAXCURV;

	// Init conditions
	int proc = 0;
	int iter = 0;
	int fail = 1;

	// Calculate score of initial curve
	cur = *c;
	calcCurveParams(e, w, h, &cur, FC_STEP);
	best = cur;

	while (proc < 2 && iter < maxIter) {
		// Use max as current
		max = cur;
		n = cur;

		if (proc == 0) {
			// Find best neighbor for (x1,y1)
			for (n.y1 = cur.y1 - FC_XYRANGE; n.y1 <= cur.y1 + FC_XYRANGE;
					n.y1 += FC_XYSTEP)	// for all neighboring y1
							{
				if (n.y1 >= 0 && n.y1 < h)	// if inside range
						{
					dy = n.y2 - n.y1;
					dy = dy * dy;

					for (n.x1 = cur.x1 - FC_XYRANGE;
							n.x1 <= cur.x1 + FC_XYRANGE; n.x1 += FC_XYSTEP)	// for all neighboring x1
									{
						if (n.x1 >= 0 && n.x1 < w)	// if inside range
								{
							dx = n.x2 - n.x1;
							dx = dx * dx;

							// If not masked
							i = n.y1 * w + n.x1;
							if (!(mask && mask[i])) {
								for (n.h = cur.h - FC_HRANGE;
										n.h <= cur.h + FC_HRANGE; n.h +=
										FC_HSTEP) // for all neighboring h
												{
									// If curvature condition is satisfied
									if (n.h * n.h < maxcurv * (dx + dy)) {
										// Calculate score
										calcCurveParams(e, w, h, &n, FC_STEP);
										if (n.score > max.score) {
											max = n;
										}
									}
								}
							}
						}
					}
				}
			}
		} else if (proc == 1) {
			// Find best neighbor for (x2,y2)
			for (n.y2 = cur.y2 - FC_XYRANGE; n.y2 <= cur.y2 + FC_XYRANGE;
					n.y2 += FC_XYSTEP)	// for all neighboring y2
							{
				if (n.y2 >= 0 && n.y2 < h)	// if inside range
						{
					dy = n.y2 - n.y1;
					dy = dy * dy;

					for (n.x2 = cur.x2 - FC_XYRANGE;
							n.x2 <= cur.x2 + FC_XYRANGE; n.x2 += FC_XYSTEP)	// for all neighboring x2
									{
						if (n.x2 >= 0 && n.x2 < w)	// if inside range
								{
							dx = n.x2 - n.x1;
							dx = dx * dx;

							// If not masked
							i = n.y1 * w + n.x1;
							if (!(mask && mask[i])) {
								for (n.h = cur.h - FC_HRANGE;
										n.h <= cur.h + FC_HRANGE; n.h +=
										FC_HSTEP) // for all neighboring h
												{
									// If curvature condition is satisfied
									if (n.h * n.h < maxcurv * (dx + dy)) {
										// Calculate score
										calcCurveParams(e, w, h, &n, FC_STEP);
										if (n.score > max.score) {
											max = n;
										}
									}
								}
							}
						}
					}
				}
			}
		}
		cur = max;

		// Take current curve if better
		if (best.score < cur.score) {
			best = cur;
			fail = 0;
		} else {
			// Go on next step
			proc++;
		}

		// Increment iteration
		iter++;
	}

	if (!fail) {
		*c = best;
	}

	return fail;
}

int fitCurve2(localEnergy *e, const int w, const int h, curve *c,
		const int maxIter, char *mask) {
	int i, dx, dy;
	curve best, cur, max, n;

	const double maxcurv = FC_MAXCURV * FC_MAXCURV;

	// Init conditions
	int cont = 1;
	int iter = 0;
	int fail = 1;

	// Calculate score of initial curve
	cur = *c;
	calcCurveParams(e, w, h, &cur, FC_STEP);
	best = cur;

	while (cont && iter < maxIter) {
		// Use max as current
		max = cur;

		// Find best neighbor for (x1,y1)
		n = cur;
		for (n.y1 = cur.y1 - FC_XYRANGE; n.y1 <= cur.y1 + FC_XYRANGE; n.y1 +=
		FC_XYSTEP)	// for all neighboring y1
				{
			if (n.y1 >= 0 && n.y1 < h)	// if inside range
					{
				for (n.y2 = cur.y2 - FC_XYRANGE; n.y2 <= cur.y2 + FC_XYRANGE;
						n.y2 += FC_XYSTEP)	// for all neighboring y2
								{
					if (n.y2 >= 0 && n.y2 < h)	// if inside range
							{
						dy = n.y2 - n.y1;
						dy = dy * dy;

						for (n.x1 = cur.x1 - FC_XYRANGE;
								n.x1 <= cur.x1 + FC_XYRANGE; n.x1 += FC_XYSTEP)	// for all neighboring x1
										{
							if (n.x1 >= 0 && n.x1 < w)	// if inside range
									{
								for (n.x2 = cur.x2 - FC_XYRANGE;
										n.x2 <= cur.x2 + FC_XYRANGE; n.x2 +=
										FC_XYSTEP)	// for all neighboring x2
												{
									if (n.x2 >= 0 && n.x2 < w)// if inside range
											{
										dx = n.x2 - n.x1;
										dx = dx * dx;

										// If not masked
										i = n.y1 * w + n.x1;
										if (!(mask && mask[i])) {
											for (n.h = cur.h - FC_HRANGE;
													n.h <= cur.h + FC_HRANGE;
													n.h += FC_HSTEP) // for all neighboring h
															{
												// If curvature condition is satisfied
												if (n.h * n.h
														< maxcurv * (dx + dy)) {
													// Calculate score
													calcCurveParams(e, w, h, &n,
													FC_STEP);
													if (n.score > max.score) {
														max = n;
													}
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
		cur = max;

		// Take current curve if better
		cont = 0;
		if (best.score < cur.score) {
			best = cur;
			cont = 1;
			fail = 0;
		}

		// Increment iteration
		iter++;
	}

	if (!fail) {
		*c = best;
	}

	return fail;
}

int fitCurve(localEnergy *e, const int w, const int h, curve *c,
		const int maxIter, char *mask) {
	int i, dx, dy;
	curve best, max, cur, n;
	float sc;
	const double maxcurv = FC_MAXCURV * FC_MAXCURV;

	// Init visited table
	int vsize = sizeof(char) * w * h;
	char *visited = (char *) malloc(vsize);
	memset(visited, 0, vsize);

	// Init conditions
	int iter = 0;
	int fail = 1;
	int cont = 1;
	int progress;

	// Init best fit curve
	best.score = 0;

	// Visit current node
	cur = *c;
	n = cur;
	i = cur.y2 * w + cur.x2;
	visited[i] = (char) 255;

	// Main loop
	while (cont && iter < maxIter) {
		// Find best neighbor
		max = cur;
		sc = -99999999.9f;
		progress = 0;
		for (n.y2 = cur.y2 - 1; n.y2 <= cur.y2 + 1; n.y2++) {
			if (n.y2 >= 0 && n.y2 < h) {
				dy = n.y2 - n.y1;
				dy = dy * dy;

				for (n.x2 = cur.x2 - 1; n.x2 <= cur.x2 + 1; n.x2++) {
					if (n.x2 >= 0 && n.x2 < w) {
						dx = n.x2 - n.x1;
						dx = dx * dx;

						// If not visited before
						i = n.y2 * w + n.x2;
						if (!visited[i] && !(mask && mask[i])) {
							for (n.h = cur.h - 1; n.h <= cur.h + 1; n.h++) {
								// If it is a new candidate
								if (cur.x2 != n.x2 || cur.y2 != n.y2
										|| cur.h != n.h) {
									// If curvature condition is satisfied
									if (n.h * n.h < maxcurv * (dx + dy)) {
										// Calculate score
										calcCurveParams(e, w, h, &n, FC_STEP);
										if (n.score > sc) {
											sc = n.score;
											max = n;
											progress = 1;
										}
									}
								}
							}
						}
					}
				}
			}
		}

		// Update visited table
		if (max.x2 != cur.x2 || max.y2 != cur.y2) {
			i = max.y2 * w + max.x2;
			visited[i] = (char) 255;
		}

		// Move to the new canditate
		cur = max;

		// Update best curve state
		if (best.score <= cur.score) {
			fail = 0;
			best = cur;
		}

		// Continuing condition
		cont = 1;
		if (best.score > EPS) {
			cont = cont && ((best.score - cur.score) / best.score <= FC_TOL); // Tolerating condition satisfied
		}
		cont = cont && progress; // and new candidate found

		// Increment iteration
		iter++;
	}

	// Return best fit curve
	if (!fail) {
		*c = best;
	}

	return fail;
}

float energystr(localEnergy *e) {
	const int r = (LE_BINS / 4) - 1;
	float se = 0;
	int j;
	for (j = e->major - r; j <= e->major + r; j++) {
		se += e->energy[j];
	}
	return se;
}

int fitCurves(localEnergy *e, const int w, const int h, curve **clist,
		const int FCS_MAXCURVE, const int FCS_XYSTEP, const int FCS_XYRANGE,
		const float FCS_INIT_THRESH, const int FCS_MINLEN) {

	int nc = 0;
	int x, y, i, j, xi, yi, xm, ym, x1, x2, y1, y2;
	float sm, si, gm;
	curve c;
	char *scanned = (char *) malloc(sizeof(char) * w * h);
	memset(scanned, 0, sizeof(char) * w * h);

	// Free previous list if exists
	if (*clist) {
		free(*clist);
	}

	// Alloc mem. for curves
	*clist = (curve *) malloc(sizeof(curve) * FCS_MAXCURVE);

	// Find global max
	gm = 0;
	for (i = 0; i < w * h; i++) {
		si = energystr(&e[i]);
		if (si > gm)
			gm = si;
	}

	// Scan blocks for local max.
	for (y = 0; y < h; y += FCS_XYSTEP) {
		y1 = y + FCS_XYSTEP / 2 - FCS_XYRANGE;
		if (y1 < 0)
			y1 = 0;
		y2 = y + FCS_XYSTEP / 2 + FCS_XYRANGE;
		if (y2 > h)
			y2 = h;

		for (x = 0; x < w; x += FCS_XYSTEP) {
			x1 = x + FCS_XYSTEP / 2 - FCS_XYRANGE;
			if (x1 < 0)
				x1 = 0;
			x2 = x + FCS_XYSTEP / 2 + FCS_XYRANGE;
			if (x2 > w)
				x2 = w;

			xm = x1;
			ym = y1;
			i = ym * w + xm;
			sm = energystr(&e[i]);
			for (yi = y1; yi < y2; yi++) {
				for (xi = x1; xi < x2; xi++) {
					i = yi * w + xi;
					si = energystr(&e[i]);
					if (si > sm) {
						sm = si;
						xm = xi;
						ym = yi;
					}
				}
			}

			i = ym * w + xm;
			if (!scanned[i] && sm / gm > FCS_INIT_THRESH) // If (not visited) and strong enough
					{
				scanned[i] = (char) 255;
				c.x1 = c.x2 = xm;
				c.y1 = c.y2 = ym;
				c.h = 0;
				if (!fitCurve3(e, w, h, &c, FC_MAXITER)) {
					if (c.len >= FCS_MINLEN)	// If curve is long enough
							{
						(*clist)[nc] = c;
						nc++;
					}
				}
			}
		}
	}

	// Clean up curve array for duplicates
	for (i = 0; i < nc - 1; i++) {
		for (j = i + 1; j < nc; j++) {
			if ((*clist)[i].x1 == (*clist)[j].x1
					&& (*clist)[i].x2 == (*clist)[j].x2
					&& (*clist)[i].y1 == (*clist)[j].y1
					&& (*clist)[i].y2 == (*clist)[j].y2
					&& (*clist)[i].h == (*clist)[j].h) {
				(*clist)[j] = (*clist)[nc - 1];
				nc--;
				j--;
			}
		}
	}

	// The other stuff
	for (i = 0; i < nc; i++) {
		// Correct sign of h param (make positive)
		if ((*clist)[i].h < 0) {
			(*clist)[i].h *= -1;
			(*clist)[i].a *= -1;
			(*clist)[i].b *= -1;
			(*clist)[i].cost *= -1;
			(*clist)[i].sint *= -1;
			x = (*clist)[i].x1;
			y = (*clist)[i].y1;
			(*clist)[i].x1 = (*clist)[i].x2;
			(*clist)[i].y1 = (*clist)[i].y2;
			(*clist)[i].x2 = x;
			(*clist)[i].y2 = y;
		}

		// Calculate aux. curve params
		calcAuxCurveParams(&((*clist)[i]));
	}

	free(scanned);
	return nc;
}

void visualizeCurves(IplImage *bckgnd, curve *clist, const int nc,
		const int ssize, IplImage *visual, const float FLC_THRESH,
		const float FLC_AVGTHRESH) {
	int i, s;
	CvScalar cl;
	cvMerge(bckgnd, bckgnd, bckgnd, NULL, visual);
	for (i = 0; i < nc; i++) {
		if (clist[i].score > FLC_THRESH && clist[i].avgscore > FLC_AVGTHRESH) {
			cl = cvScalar(255, 0, 0);
			s = 2;
		} else {
			cl = cvScalar(0, 0, 255);
			s = 1;
		}
		drawCurve1(visual, &clist[i], ssize, s, FC_STEP, cl);
	}
}

void initCurveFitting() {
	int j;
	double t;
	for (j = 0; j < LE_BINS; j++) {
		t = (PI * j) / LE_BINS;
		angle_x[j] = cos(t);
		angle_y[j] = sin(t);
	}
}

void drawBinaryCurves(curve *clist, const int nc, const int psize,
		IplImage **out, const int w, const int h) {
	int i, r = psize / 2;

	if (*out) {
		cvReleaseImage(out);
	}

	*out = cvCreateImage(cvSize(w, h), 8, 1);
	cvSetZero(*out);

	for (i = 0; i < nc; i++) {
		drawCurve(*out, &clist[i], 1, r, FC_STEP, cvScalar(255));
	}
}

void visualizeHiLoCurves(curve *clist_lo, const int nc_lo, const int w_lo,
		const int h_lo, curve *clist_hi, const int nc_hi,
		const double maxoverlap, IplImage *bckgnd, IplImage *visual) {
	IplImage *memtemp = 0;
	// Prepare low freq. curves coverage
	drawBinaryCurves(clist_lo, nc_lo, (int) CURV_HILO_COV, &memtemp, w_lo,
			h_lo);

	cvMerge(bckgnd, bckgnd, bckgnd, NULL, visual);
	double scale = (double) LE_SSIZE_LO / LE_SSIZE_HI;
	int i, xi, yi, in, tot;
	curve *c;
	double t, k, x, y, m, dt;

	for (i = 0; i < nc_hi; i++) {
		c = &clist_hi[i];

		in = 0;
		tot = 0;

		t = 0;
		while (t <= c->R && c->R >= EPS) {
			k = c->a * t * t + c->b * t;
			x = c->cost * t - c->sint * k + c->x1;
			y = c->sint * t + c->cost * k + c->y1;

			xi = cvRound(x * scale);
			yi = cvRound(y * scale);
			if (xi >= 0 && yi >= 0 && xi < w_lo && yi < h_lo) {
				if (cvGetReal2D(memtemp, yi, xi))
					in++;
			}
			tot++;

			m = 2 * c->a * t + c->b;
			dt = FC_STEP / sqrt(1.0 + m * m);
			t += dt;
		}

		if ((double) in / tot < maxoverlap) {
			drawCurve1(visual, c, (int) LE_SSIZE_HI, 2, FC_STEP,
					cvScalar(0, 0, 255));
		}

	}

	for (i = 0; i < nc_lo; i++) {
		drawCurve1(visual, &clist_lo[i], (int) LE_SSIZE_LO, 2, FC_STEP,
				cvScalar(255, 0, 0));
	}

	cvReleaseImage(&memtemp);
}

void getHixCurves(curve *clist_lo, const int nc_lo, const int w_lo,
		const int h_lo, curve *clist_hi, const int nc_hi, curve **clist_hix,
		int *nc_hix) {
	*clist_hix = (curve *) malloc(sizeof(curve) * nc_hi);
	*nc_hix = 0;

	IplImage *memtemp = 0;
	// Prepare low freq. curves coverage
	drawBinaryCurves(clist_lo, nc_lo, (int) CURV_HILO_COV, &memtemp, w_lo,
			h_lo);

	double scale = (double) LE_SSIZE_LO / LE_SSIZE_HI;
	int i, xi, yi, in, tot;
	curve *c;
	double t, k, x, y, m, dt;

	for (i = 0; i < nc_hi; i++) {
		c = &clist_hi[i];

		in = 0;
		tot = 0;

		t = 0;
		while (t <= c->R && c->R >= EPS) {
			k = c->a * t * t + c->b * t;
			x = c->cost * t - c->sint * k + c->x1;
			y = c->sint * t + c->cost * k + c->y1;

			xi = cvRound(x * scale);
			yi = cvRound(y * scale);
			if (xi >= 0 && yi >= 0 && xi < w_lo && yi < h_lo) {
				if (cvGetReal2D(memtemp, yi, xi))
					in++;
			}
			tot++;

			m = 2 * c->a * t + c->b;
			dt = FC_STEP / sqrt(1.0 + m * m);
			t += dt;
		}

		if ((double) in / tot < CURV_HILO_COV_TH) {
			(*clist_hix)[*nc_hix] = clist_hi[i];
			(*nc_hix)++;
		}
	}

	cvReleaseImage(&memtemp);
}

void filterCurves(const curve *cin, const int ncin, curve **cout, int *ncout,
		const float FLC_THRESH, const float FLC_AVGTHRESH) {
	int i;
	*cout = (curve *) malloc(sizeof(curve) * ncin);
	*ncout = 0;
	for (i = 0; i < ncin; i++) {
		if (cin[i].score > FLC_THRESH && cin[i].avgscore > FLC_AVGTHRESH) {
			(*cout)[*ncout] = cin[i];
			(*ncout)++;
		}
	}
}

void getCurvePoints(const curve *c, const double S, double *cx, double *cy) {
	int j;
	double t, k, dt, m, x, y;

	// Init point counter
	j = 0;
	// Scan "len" points (0<=t<=R)
	t = 0;
	while (j < c->len) {
		k = c->a * t * t + c->b * t;
		x = c->cost * t - c->sint * k + c->x1;
		y = c->sint * t + c->cost * k + c->y1;

		cx[j] = x;
		cy[j] = y;

		m = 2 * c->a * t + c->b;
		dt = S / sqrt(1.0 + m * m);
		t += dt;
		j++;
	}
}

void saveLocalEnergyStructure(int s, char const *tag, localEnergy *e, int w,
		int h) {
	FILE *f;
	char outputfn[256];
	char outputsrcfn[256];
	sprintf(outputsrcfn, "%s%s%s.dat", DESTPATH, tag, FNAME);
	sprintf(outputfn, outputsrcfn, s);
	if (e) {
		printf("Saving: %s\n", outputfn);
		f = fopen(outputfn, "wb");
		fwrite(&w, sizeof(int), 1, f);
		fwrite(&h, sizeof(int), 1, f);
		fwrite(e, sizeof(localEnergy), w * h, f);
		fclose(f);
	}
}

void loadLocalEnergyStructure(int s, char const *tag, localEnergy **e, int *w,
		int *h) {
	FILE *f;
	size_t t;
	char outputfn[256];
	char outputsrcfn[256];
	sprintf(outputsrcfn, "%s%s%s.dat", DESTPATH, tag, FNAME);
	sprintf(outputfn, outputsrcfn, s);
	if (*e) {
		free(*e);
	}
	printf("Loading: %s\n", outputfn);
	f = fopen(outputfn, "rb");
	t = fread(w, sizeof(int), 1, f);
	t = fread(h, sizeof(int), 1, f);
	*e = (localEnergy *) malloc(sizeof(localEnergy) * (*w) * (*h));
	t = fread(*e, sizeof(localEnergy), (*w) * (*h), f);
	if (!t)
		printf("Error loading data!\n");
	fclose(f);
}

void saveCurveStructure(int s, char const *tag, curve *c, int numRecords) {
	FILE *f;
	char outputfn[256];
	char outputsrcfn[256];
	sprintf(outputsrcfn, "%s%s%s.dat", DESTPATH, tag, FNAME);
	sprintf(outputfn, outputsrcfn, s);
	if (c) {
		printf("Saving: %s\n", outputfn);
		f = fopen(outputfn, "wb");
		fwrite(&numRecords, sizeof(int), 1, f);
		fwrite(c, sizeof(curve), numRecords, f);
		fclose(f);
	}
}

int loadCurveStructure(int s, char const *tag, curve **c) {
	int numRecords;
	FILE *f;
	size_t t;
	char outputfn[256];
	char outputsrcfn[256];
	sprintf(outputsrcfn, "%s%s%s.dat", DESTPATH, tag, FNAME);
	sprintf(outputfn, outputsrcfn, s);
	if (*c) {
		free(*c);
	}
	printf("Loading: %s\n", outputfn);
	f = fopen(outputfn, "rb");
	t = fread(&numRecords, sizeof(int), 1, f);
	*c = (curve *) malloc(sizeof(curve) * numRecords);
	t = fread(*c, sizeof(curve), numRecords, f);
	fclose(f);
	if (!t)
		printf("Error loading data!\n");
	return numRecords;
}

int isInside(const curve *c, double x, double y, double bias) {
	x -= c->x1;
	y -= c->y1;
	double t = c->cost * x + c->sint * y;
	double k = -c->sint * x + c->cost * y;
	if (k + bias <= c->a * t * t + c->b * t)
		return 1;
	return 0;
}
