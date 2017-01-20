/*
 * poly.cpp
 *
 *  Created on: Jan 15, 2017
 *      Author: fst
 */

#include "poly.h"

void initBlankPoly2d(int n, poly2d *p) {
	if (p) {
		p->next = NULL;
		p->n = n;
		if (n) {
			p->vertex = (CvPoint *) malloc(sizeof(CvPoint) * n);
			memset(p->vertex, 0, sizeof(CvPoint) * n);
		} else
			p->vertex = NULL;
	}
}

void initPoly2d(double (*v)[2], int n, poly2d *p) {
	if (p) {
		poly2d t;

		initBlankPoly2d(n, &t);
		for (int i = 0; i < n; i++) {
			t.vertex[i].x = (int) v[i][0];
			t.vertex[i].y = (int) v[i][1];
		}
		if (p->next) {
			deinitPoly2d(p->next);
			free(p->next);
		}
		p->next = NULL;
		combinePoly2d(&t, NULL, 0, p);
	}
}

void deinitPoly2d(poly2d *p) {
	if (p) {
		if (p->vertex)
			free(p->vertex);
		p->vertex = NULL;
		p->n = 0;
		if (p->next) {
			deinitPoly2d(p->next);
			free(p->next);
		}
		p->next = NULL;
	}
}

void initBlankPoly25d(int t, int start_z, poly25d *p) {
	if (p) {
		p->t = t;
		p->start_z = start_z;
		p->end_z = start_z + t - 1;
		if (t) {
			p->slice = (poly2d *) malloc(sizeof(poly2d) * t);
			memset(p->slice, 0, sizeof(poly2d) * t);
		} else
			p->slice = NULL;
	}
}

void deinitPoly25d(poly25d *p) {
	if (p) {
		if (p->slice) {
			for (int i = 0; i < p->t; i++)
				deinitPoly2d(&p->slice[i]);
			free(p->slice);
			p->slice = NULL;
		}
		p->t = 0;
		p->start_z = p->end_z = -1;
	}
}

void initPoly25d(snake25d *s, poly25d *p) {
	if (s && p) {
		initBlankPoly25d(s->t, s->start_z, p);
		for (int i = 0; i < s->t; i++) {
			initPoly2d(s->node[i], SN_N, &p->slice[i]);
		}
	}
}

void initPoly25dArray(snake25d *sarray, int n, poly25d *parray) {
	if (sarray && parray) {
		for (int i = 0; i < n; i++) {
			initPoly25d(&sarray[i], &parray[i]);
		}
	}
}

void deinitPoly25dArray(poly25d *parray, int n) {
	if (parray) {
		for (int i = 0; i < n; i++) {
			deinitPoly25d(&parray[i]);
		}
	}
}

int poly2dArea(poly2d *p, int *r) {
	int k, j;
	int area = 0;
	for (k = 0, j = 1; k < p->n; k++, j = (j + 1) % p->n) {
		area += p->vertex[k].x * p->vertex[j].y
				- p->vertex[k].y * p->vertex[j].x;
	}
	area /= 2;
	if (area < 0) {
		if (r)
			*r = 1;
		area = -area;
	} else {
		if (r)
			*r = 0;
	}
	return area;
}

void correctPoly2d(poly2d *p) {
	int i, j, r = 0;
	CvPoint t;
	poly2dArea(p, &r);
	if (r) {
		for (i = 0, j = p->n - 1; i < p->n / 2; i++, j--) {
			t = p->vertex[i];
			p->vertex[i] = p->vertex[j];
			p->vertex[j] = t;
		}
	}
	if (p->next)
		correctPoly2d(p->next);
}

void combinePoly2d(poly2d *p1, poly2d *p2, int intersect, poly2d *out, int *a1,
		int *a2, int *a12) {
	int w = 0, h = 0;
	int i;
	CvPoint *k1;
	CvPoint *k2;
	int n1, n2;
	poly2d *p1_seq;
	poly2d *p2_seq;

	// Determine working area
	p1_seq = p1;
	p2_seq = p2;
	do {
		if (p1_seq && p1_seq->n) {
			n1 = p1_seq->n;
			k1 = p1_seq->vertex;
		} else {
			n1 = 0;
			k1 = NULL;
		}
		if (p2_seq && p2_seq->n) {
			n2 = p2_seq->n;
			k2 = p2_seq->vertex;
		} else {
			n2 = 0;
			k2 = NULL;
		}

		for (i = 0; i < n1; i++) {
			if (w < k1[i].x)
				w = k1[i].x;
			if (h < k1[i].y)
				h = k1[i].y;
		}
		for (i = 0; i < n2; i++) {
			if (w < k2[i].x)
				w = k2[i].x;
			if (h < k2[i].y)
				h = k2[i].y;
		}

		// Go to next poly
		if (p1_seq)
			p1_seq = p1_seq->next;
		if (p2_seq)
			p2_seq = p2_seq->next;
	} while (p1_seq || p2_seq);
	w += 2;
	h += 2;

	// Initialize buffers
	IplImage *buf1 = cvCreateImage(cvSize(w, h), IPL_DEPTH_8U, 1);
	IplImage *buf2 = cvCreateImage(cvSize(w, h), IPL_DEPTH_8U, 1);
	cvSetZero(buf1);
	cvSetZero(buf2);

	// Draw polygons
	p1_seq = p1;
	p2_seq = p2;
	do {
		if (p1_seq && p1_seq->n) {
			n1 = p1_seq->n;
			k1 = p1_seq->vertex;
		} else {
			n1 = 0;
			k1 = NULL;
		}
		if (p2_seq && p2_seq->n) {
			n2 = p2_seq->n;
			k2 = p2_seq->vertex;
		} else {
			n2 = 0;
			k2 = NULL;
		}

		if (n1)
			cvFillPoly(buf1, &k1, &n1, 1, cvScalar(255, 255, 255));
		if (n2)
			cvFillPoly(buf2, &k2, &n2, 1, cvScalar(255, 255, 255));

		// Go to next poly
		if (p1_seq)
			p1_seq = p1_seq->next;
		if (p2_seq)
			p2_seq = p2_seq->next;
	} while (p1_seq || p2_seq);

	// Obtain area of the polygons
	if (a1)
		*a1 = cvCountNonZero(buf1);
	if (a2)
		*a2 = cvCountNonZero(buf2);

	// Create combined polygon
	if (intersect)
		cvAnd(buf1, buf2, buf1);
	else
		cvOr(buf1, buf2, buf1);

	// Obtain combined polygon area
	if (a12)
		*a12 = cvCountNonZero(buf1);

	// Extract contours
	if (out) {
		CvMemStorage *storage = cvCreateMemStorage(0);
		CvSeq *contour = NULL;
		cvFindContours(buf1, storage, &contour, sizeof(CvContour),
				CV_RETR_EXTERNAL, CV_CHAIN_APPROX_NONE);

		// Prepare contour output
		CvSeq *seq = contour;
		poly2d *polySeq = out, *prev = NULL;
		deinitPoly2d(out);
		while (seq) {
			if (seq->first->count) {
				initBlankPoly2d(seq->first->count, polySeq);
				memcpy(polySeq->vertex, seq->first->data,
						sizeof(CvPoint) * polySeq->n);
				polySeq->next = NULL;
				prev = polySeq;
				polySeq = NULL;
			}
			seq = seq->h_next;
			if (seq && seq->first->count) {
				if (!polySeq) {
					polySeq = prev->next = (poly2d *) calloc(1, sizeof(poly2d));
					deinitPoly2d(polySeq);
				}
			}
		}
		correctPoly2d(out);

		cvReleaseMemStorage(&storage);
	}

	// Free memory
	cvReleaseImage(&buf1);
	cvReleaseImage(&buf2);
}

void combinePoly25d(poly25d *p1, poly25d *p2, int intersect, poly25d *out,
		int *a1, int *a2, int *a12) {
	int start_z =
			(((p1->start_z < p2->start_z) & 1) ^ intersect) ?
					p1->start_z : p2->start_z;
	int end_z =
			(((p1->end_z > p2->end_z) & 1) ^ intersect) ? p1->end_z : p2->end_z;

	poly2d *s1, *s2, *pout;
	int t = end_z - start_z + 1;
	int z;
	int b1, b2, b12, c1 = 0, c2 = 0, c12 = 0;

	if (out)
		deinitPoly25d(out);
	if (t > 0) {
		if (out)
			initBlankPoly25d(t, start_z, out);

		for (z = start_z; z <= end_z; z++) {
			if (p1->start_z <= z && p1->end_z >= z)
				s1 = &p1->slice[z - p1->start_z];
			else
				s1 = NULL;
			if (p2->start_z <= z && p2->end_z >= z)
				s2 = &p2->slice[z - p2->start_z];
			else
				s2 = NULL;

			if (out)
				pout = &out->slice[z - start_z];
			else
				pout = NULL;

			combinePoly2d(s1, s2, intersect, pout, &b1, &b2, &b12);
			c1 += b1;
			c2 += b2;
			c12 += b12;
		}
	}

	if (a1)
		*a1 = c1;
	if (a2)
		*a2 = c2;
	if (a12)
		*a12 = c12;
}

void safecopyPoly2d(poly2d *source, poly2d *dest) {
	if (source && dest) {
		deinitPoly2d(dest);
		initBlankPoly2d(source->n, dest);
		memcpy(dest->vertex, source->vertex, sizeof(CvPoint) * source->n);
		dest->next = NULL;
		if (source->next) {
			dest->next = (poly2d *) calloc(1, sizeof(poly2d));
			safecopyPoly2d(source->next, dest->next);
		}
	}
}

void safecopyPoly25d(poly25d *source, poly25d *dest) {
	if (source && dest) {
		deinitPoly25d(dest);
		initBlankPoly25d(source->t, source->start_z, dest);
		for (int i = 0; i < source->t; i++) {
			safecopyPoly2d(&source->slice[i], &dest->slice[i]);
		}
	}
}

int mergeArrayOfPoly25d(poly25d *parray, int n, double th, poly25d *outArray) {
	int n_out = 0;
	int a1, a2, a12;
	int i, j;
	int cont_merge = 1;
	double d1, d2;

	poly25d *temp = (poly25d *) calloc(1, sizeof(poly25d));
	deinitPoly25dArray(outArray, n);

	if (parray && n) {
		for (i = 0; i < n; i++)
			safecopyPoly25d(&parray[i], &outArray[i]);
		n_out = n;

		while (cont_merge) {
			cont_merge = 0;

			for (i = 0; i < n_out - 1; i++) {
				printf("\rMerging: %d / %d (thres = %f)          ", i,
						n_out - 1, th);
				fflush(stdout);
				for (j = i + 1; j < n_out; j++) {
					a1 = a2 = a12 = 0;
					combinePoly25d(&outArray[i], &outArray[j], 1, NULL, &a1,
							&a2, &a12);
					if (a1)
						d1 = (double) a12 / a1;
					else
						d1 = 0;
					if (a2)
						d2 = (double) a12 / a2;
					else
						d2 = 0;
					if (d1 >= th || d2 >= th) // Merging condition
							{
						combinePoly25d(&outArray[i], &outArray[j], 0, temp);
						safecopyPoly25d(&outArray[n_out - 1], &outArray[j]);
						safecopyPoly25d(temp, &outArray[i]);
						deinitPoly25d(&outArray[n_out - 1]);
						n_out--;
						cont_merge = 1;
					}
				}
			}
		}
		printf("\n");
	}

	// Free mem.
	deinitPoly25d(temp);
	free(temp);

	return n_out;
}

int convertValidSnakesToPolyArray(snake25d *sarray, int n, double th,
		poly25d **parray) {
	snake25d *vs = NULL;
	int nvs = filterSnakeArrayByValidity(sarray, n, th, &vs);
	*parray = (poly25d *) calloc(nvs, sizeof(poly25d));
	initPoly25dArray(vs, nvs, *parray);
	if (vs)
		free(vs);
	return nvs;
}

int getTriangles(poly2d *p1, poly2d *p2, int (**faces)[3]) {
	if (!p1 || !p2)
		return 0;

	if (*faces)
		free(*faces);
	*faces = (int (*)[3]) calloc(p1->n + p2->n, sizeof(int[3]));
	int n_faces = 0;

	int i, j, ii, jj, k;
	int dx, dy, d;
	int dmin = 9999999;
	int is = 0, js = 0;

	for (i = 0; i < p1->n; i++) {
		for (j = 0; j < p2->n; j++) {
			d = 0;
			ii = i;
			jj = j;
			for (k = 0; k < 20; k++) {
				dx = p1->vertex[ii].x - p2->vertex[jj].x;
				dy = p1->vertex[ii].y - p2->vertex[jj].y;
				d += 1 + dx * dx + dy * dy;
				ii++;
				if (ii >= p1->n)
					ii -= p1->n;
				jj++;
				if (jj >= p2->n)
					jj -= p2->n;
			}
			if (d < dmin) {
				dmin = d;
				is = i + 10;
				if (is >= p1->n)
					is -= p1->n;
				js = j + 10;
				if (js >= p2->n)
					js -= p2->n;
			}
		}
	}

	double id, jd;
	i = j = 0;
	while (n_faces < p1->n + p2->n) {
		id = (double) i / p1->n;
		jd = (double) j / p2->n;
		ii = (i + is) % p1->n;
		jj = (j + js) % p2->n;
		if (id < jd) {
			(*faces)[n_faces][0] = ii;
			(*faces)[n_faces][1] = (ii + 1) % p1->n;
			(*faces)[n_faces][2] = jj + p1->n;
			i++;
		} else {
			(*faces)[n_faces][0] = ii;
			(*faces)[n_faces][2] = jj + p1->n;
			(*faces)[n_faces][1] = (jj + 1) % p2->n + p1->n;
			j++;
		}
		n_faces++;
	}

	return n_faces;
}

void savePolyArrayAsPLY(poly25d *p, int n) {
	double xy_scale = LE_SSIZE_LO * TFACTOR;
	double z_scale = RESOLUTION;
	double x, y, z;
	poly2d *seq;
	int k, i, j, n_vertex, cv;
	FILE *f;
	int (*faces)[3] = NULL;
	int n_faces;

	if (p) {
		n_vertex = 0;
		for (k = 0; k < n; k++) {
			for (i = 0; i < p[k].t; i++) {
				seq = &p[k].slice[i];
				while (seq) {
					n_vertex += seq->n;
					seq = seq->next;
				}
			}
		}

		n_faces = 0;
		for (k = 0; k < n; k++) {
			for (i = 0; i < p[k].t; i++) {
				if (i < p[k].t - 1) {
					n_faces += getTriangles(&p[k].slice[i], &p[k].slice[i + 1],
							&faces);
				}
			}
		}

		char outputfn[256];
		char outputsrcfn[256];
		sprintf(outputsrcfn, "%s%s%s.ply", DESTPATH, "poly_", FNAME);
		sprintf(outputfn, outputsrcfn, SLICE_START);
		printf("Saving: %s\n", outputfn);
		f = fopen(outputfn, "w");

		fprintf(f, "ply\nformat ascii 1.0\n");
		fprintf(f,
				"element vertex %d\nproperty float x\nproperty float y\nproperty float z\n",
				n_vertex);
		fprintf(f, "element face %d\nproperty list uchar int vertex_indices\n",
				n_faces);
		fprintf(f, "end_header\n");

		for (k = 0; k < n; k++) {
			for (i = 0; i < p[k].t; i++) {
				z = (i + p[k].start_z) * z_scale;
				seq = &p[k].slice[i];
				while (seq) {
					for (int j = 0; j < seq->n; j++) {
						x = seq->vertex[j].x * xy_scale;
						y = seq->vertex[j].y * xy_scale;
						fprintf(f, "%f %f %f\n", x, y, z);
					}
					seq = seq->next;
				}
			}
		}

		cv = 0;
		for (k = 0; k < n; k++) {
			for (i = 0; i < p[k].t; i++) {
				seq = &p[k].slice[i];
				if (i < p[k].t - 1) {
					n_faces = getTriangles(&p[k].slice[i], &p[k].slice[i + 1],
							&faces);
					for (j = 0; j < n_faces; j++) {
						fprintf(f, "3 %d %d %d\n", faces[j][0] + cv,
								faces[j][1] + cv, faces[j][2] + cv);
					}
				}
				while (seq) {
					cv += seq->n;
					seq = seq->next;
					if (seq)
						printf("Warning: Poly2d list is not supported!\n");
				}
			}
		}

		fclose(f);
	}

	if (faces)
		free(faces);
}
