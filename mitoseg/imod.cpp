/*
 * imod.cpp
 *
 *  Created on: Sep 28, 2017
 *      Author: fst
 */

#include "imod.h"

int imod_endian_reverse(int x) {
	return ((x & 0xff000000) >> 24) | ((x & 0x00ff0000) >> 8)
			| ((x & 0x0000ff00) << 8) | ((x & 0x000000ff) << 24);
}

float imod_endian_reverse(float x) {
	void *p = &x;
	int y = imod_endian_reverse(*((int *) p));
	p = &y;
	return *((float *) p);
}

void savePolyArrayAsIMOD(poly25d *p, int n, int x_size, int y_size, int z_size,
		double xy_scale, double z_scale, int x_shift, int y_shift) {
	FILE *f;
	int id, k, i, j;
	imod_model model;
	imod_objt objt;
	imod_cont cont;
	float pt;
	char outputfn[256];
	char outputsrcfn[256];
	sprintf(outputsrcfn, "%s%s%s.mod", DESTPATH, "imod_", FNAME);
	sprintf(outputfn, outputsrcfn, SLICE_START);
	printf("Saving: %s\n", outputfn);
	f = fopen(outputfn, "w");

	// Write header
	id = IMOD_VALUE(IMOD_FILE_ID);
	fwrite(&id, sizeof(int), 1, f);
	id = IMOD_VALUE(IMOD_VERSION_ID);
	fwrite(&id, sizeof(int), 1, f);
	// Write model info
	memset(&model, 0, sizeof(imod_model));
	model.objsize = IMOD_VALUE(n);
	model.xmax = IMOD_VALUE(x_size);
	model.ymax = IMOD_VALUE(y_size);
	model.zmax = IMOD_VALUE(z_size);
	model.blacklevel = IMOD_VALUE(0);
	model.whitelevel = IMOD_VALUE(255);
	model.xscale = IMOD_VALUE(1.f);
	model.yscale = IMOD_VALUE(1.f);
	model.zscale = IMOD_VALUE(1.f);
	fwrite(&model, sizeof(imod_model), 1, f);
	// Write objects' data
	for (k = 0; k < n; k++) {
		id = IMOD_VALUE(IMOD_OBJT_ID);
		fwrite(&id, sizeof(int), 1, f);
		memset(&objt, 0, sizeof(imod_objt));
		objt.contsize = IMOD_VALUE(p[k].t);
		objt.red = IMOD_VALUE(0.f);
		objt.green = IMOD_VALUE(0.f);
		objt.blue = IMOD_VALUE(1.f);
		objt.symsize = 1;
		sprintf(objt.name, "obj_%d", k + 1);
		fwrite(&objt, sizeof(imod_objt), 1, f);
		// Write object's contours
		for (i = 0; i < p[k].t; i++) {
			poly2d *ptr = &p[k].slice[i];
			while (ptr) {
				id = IMOD_VALUE(IMOD_CONT_ID);
				fwrite(&id, sizeof(int), 1, f);
				memset(&cont, 0, sizeof(imod_cont));
				cont.psize = IMOD_VALUE(ptr->n);
				fwrite(&cont, sizeof(imod_cont), 1, f);
				// Write contour points
				for (j = 0; j < ptr->n; j++) {
					pt = (float) (ptr->vertex[j].x * xy_scale + x_shift);
					pt = IMOD_VALUE(pt);
					fwrite(&pt, sizeof(float), 1, f);
					pt = (float) (y_size - 1 - ptr->vertex[j].y * xy_scale
							- y_shift);
					pt = IMOD_VALUE(pt);
					fwrite(&pt, sizeof(float), 1, f);
					pt = (float) ((p[k].start_z + i) * z_scale);
					pt = IMOD_VALUE(pt);
					fwrite(&pt, sizeof(float), 1, f);
				}
				ptr = ptr->next;
			}
		}
	}
	// EOF
	id = IMOD_VALUE(IMOD_IEOF_ID);
	fwrite(&id, sizeof(int), 1, f);
	//
	fclose(f);
}
