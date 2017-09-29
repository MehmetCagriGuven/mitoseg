/*
 * main.cpp
 *
 *  Created on: Jan 13, 2017
 *      Author: fst
 */

#define VERSION "1.1"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "datasets.h"
#include "settings.h"
#include "phase1_main.h"
#include "phase2_main.h"
#include "phase3_main.h"

bool zrangeSet = false;
bool psizeSet = false;
bool roiSet = false;
bool srcSet = false;
bool dstSet = false;
bool phaseSet = false;
bool fnameSet = false;
bool validSet = false;
bool thickSet = false;
int phaseNum = 1;
int thickParam = 0;

void displayHelp() {
	printf("MitoSeg v" VERSION " by FST\n");
	printf("Based on the paper:\n");
	printf(
			"Tasel, S.F., Mumcuoglu, E.U., Hassanpour, R.Z., Perkins, G., 2016.\nA validated active contour method driven by parabolic arc model\nfor detection and segmentation of mitochondria.\nJ. Struct. Biol. 194, 253–271. doi:10.1016/j.jsb.2016.03.002\n\n");
	printf("Usage:\n");
	printf("mitoseg -options <parameters> <filename pattern>\n\n");
	printf("Options:\n");
	printf("\t-zrange <start slice #> <end slice #>\tSpecify z-range\n");
	printf("\t-psize <pixel size>\t\t\tSpecify pixel size as nm/px\n");
	printf(
			"\t-roi <left> <top> <width> <height>\tSpecify region of interest\n");
	printf("\t-src <source directory>\t\t\tSpecify source directory\n");
	printf("\t-dst <destination directory>\t\tSpecify destination directory\n");
	printf(
			"\t-phase <phase #>\t\t\tApply phase <phase #> only.\n\t\t\t\t\t\t(must be 1, 2 or 3)\n");
	printf(
			"\t-valid <validity threshold>\t\tSpecify validity threshold\n\t\t\t\t\t\tbetween 0-1 (default: 0.75)\n");
	printf(
			"\t-thick <z-thickness>\t\t\tSpecify snake z-thickness\n\t\t\t\t\t\tbetween 5-500 (default: 20)\n\t\t\t\t\t\tset 'full' to use z-range\n");
	printf(
			"\t-cores <# of cpu cores>\t\t\tSet # of cpu cores to utilize\n\t\t\t\t\t\t(default: 1)\n");
	printf(
			"\n\nThe options -zrange, -psize and <filename pattern> are mandatory input.\n\n");
	printf("Example:\n");
	printf("\tmitoseg -zrange 30 100 -psize 2.0 dataset_slice%%04d.tif\n");
	printf(
			"\t\tProcess files dataset_slice0030.tif ... dataset_slice0100.tif\n");
	printf("\t\tassuming that pixel size is 2.0nm\n");
	printf("\tmitoseg -zrange 40 120 -psize 1.1 mito%%d.bmp\n");
	printf("\t\tProcess files mito40.bmp ... mito120.bmp\n");
	printf("\t\tassuming that pixel size is 1.1nm\n");
	printf("\nE-mail to fst@cankaya.edu.tr for more information.\n\n");
}

int main(int argc, char** argv) {

	if (argc < 2) {
		displayHelp();
		exit(0);
	}

	int c = 1;
	while (c < argc) {
		if (strcmp(argv[c], "-zrange") == 0) {
			if (c + 2 < argc) {
				SLICE_START = atoi(argv[++c]);
				SLICE_END = atoi(argv[++c]);
				zrangeSet = true;
				if (SLICE_START < 0 || SLICE_END < 0
						|| SLICE_START > SLICE_END) {
					printf("Error: Invalid parameter for -zrange\n");
					printf("Use mitoseg with no parameter for help\n\n");
					exit(-1);
				}
			} else {
				printf("Error: Insufficient parameter for -zrange\n");
				printf("Use mitoseg with no parameter for help\n\n");
				exit(-1);
			}
		} else if (strcmp(argv[c], "-psize") == 0) {
			if (c + 1 < argc) {
				RESOLUTION = atof(argv[++c]);
				psizeSet = true;
				if (RESOLUTION <= 0) {
					printf("Error: Invalid parameter for -psize\n");
					printf("Use mitoseg with no parameter for help\n\n");
					exit(-1);
				}
			} else {
				printf("Error: Insufficient parameter for -psize\n");
				printf("Use mitoseg with no parameter for help\n\n");
				exit(-1);
			}
		} else if (strcmp(argv[c], "-roi") == 0) {
			if (c + 4 < argc) {
				ROI_X = atoi(argv[++c]);
				ROI_Y = atoi(argv[++c]);
				ROI_W = atoi(argv[++c]);
				ROI_H = atoi(argv[++c]);
				roiSet = true;
				if (ROI_X < 0 || ROI_Y < 0 || ROI_W <= 0 || ROI_H <= 0) {
					printf("Error: Invalid parameter for -roi\n");
					printf("Use mitoseg with no parameter for help\n\n");
					exit(-1);
				}
			} else {
				printf("Error: Insufficient parameter for -roi\n");
				printf("Use mitoseg with no parameter for help\n\n");
				exit(-1);
			}
		} else if (strcmp(argv[c], "-src") == 0) {
			if (c + 1 < argc) {
				strcpy(SOURCEPATH, argv[++c]);
				srcSet = true;
			} else {
				printf("Error: Insufficient parameter for -src\n");
				printf("Use mitoseg with no parameter for help\n\n");
				exit(-1);
			}
		} else if (strcmp(argv[c], "-dst") == 0) {
			if (c + 1 < argc) {
				strcpy(DESTPATH, argv[++c]);
				dstSet = true;
			} else {
				printf("Error: Insufficient parameter for -dst\n");
				printf("Use mitoseg with no parameter for help\n\n");
				exit(-1);
			}
		} else if (strcmp(argv[c], "-phase") == 0) {
			if (c + 1 < argc) {
				phaseNum = atoi(argv[++c]);
				phaseSet = true;
				if (phaseNum < 1 || phaseNum > 3) {
					printf("Error: Phase # must be between 1 and 3.\n");
					printf("Use mitoseg with no parameter for help\n\n");
					exit(-1);
				}
			} else {
				printf("Error: Insufficient parameter for -phase\n");
				printf("Use mitoseg with no parameter for help\n\n");
				exit(-1);
			}
		} else if (strcmp(argv[c], "-valid") == 0) {
			if (c + 1 < argc) {
				POLY_VALIDITY = atof(argv[++c]);
				validSet = true;
				if (POLY_VALIDITY < 0 || POLY_VALIDITY > 1) {
					printf(
							"Error: Invalid parameter for -valid (must be between 0-1)\n");
					printf("Use mitoseg with no parameter for help\n\n");
					exit(-1);
				}
			} else {
				printf("Error: Insufficient parameter for -valid\n");
				printf("Use mitoseg with no parameter for help\n\n");
				exit(-1);
			}
		} else if (strcmp(argv[c], "-thick") == 0) {
			if (c + 1 < argc) {
				if (strcmp(argv[c + 1], "full") == 0) {
					thickSet = true;
				} else {
					sn25d_t = atoi(argv[c + 1]);
					if (sn25d_t < 5 || sn25d_t > 500) {
						printf(
								"Error: Invalid parameter for -thick (must be between 5-500)\n");
						printf("Use mitoseg with no parameter for help\n\n");
						exit(-1);
					}
				}
				c++;
			} else {
				printf("Error: Insufficient parameter for -thick\n");
				printf("Use mitoseg with no parameter for help\n\n");
				exit(-1);
			}
		} else if (strcmp(argv[c], "-cores") == 0) {
			if (c + 1 < argc) {
				numCores = atoi(argv[++c]);
				if (numCores < 1 || numCores > 256) {
					printf("Error: # of cores must be between 1 and 256.\n");
					printf("Use mitoseg with no parameter for help\n\n");
					exit(-1);
				}
			} else {
				printf("Error: Insufficient parameter for -cores\n");
				printf("Use mitoseg with no parameter for help\n\n");
				exit(-1);
			}
		} else if (argv[c][0] == '-') {
			printf("Error: Unknown option %s\n", argv[c]);
			printf("Use mitoseg with no parameter for help\n\n");
			exit(-1);
		} else {
			strcpy(FNAME, argv[c]);
			fnameSet = true;
		}
		c++;
	}
	if (!(zrangeSet && psizeSet && fnameSet)) {
		printf(
				"The options -zrange, -psize and <filename pattern> are mandatory input.\n");
		printf("Use mitoseg with no parameter for help\n\n");
		exit(-1);
	}
	if (thickSet) {
		sn25d_t = SLICE_END - SLICE_START + 1;
		if (sn25d_t < 5)
			sn25d_t = 5;
	}
	if (SLICE_END - SLICE_START + 1 < sn25d_t) {
		printf("Error: Z-range must contain at least %d slices.\n", sn25d_t);
		printf("Use mitoseg with no parameter for help\n\n");
		exit(-1);
	}
	printf("Using pixel size: %.2fnm\n", RESOLUTION);
	if (!roiSet) {
		printf("ROI was not set. Enabling Auto-ROI...\n");
		ROI_X = ROI_Y = ROI_W = ROI_H = -1;
	} else
		printf("Using ROI: %d %d %d %d\n", ROI_X, ROI_Y, ROI_W, ROI_H);
	printf("Using snake z-thickness = %d\n", sn25d_t);
	printf("Using validity threshold = %.2f\n", POLY_VALIDITY);
	if (!srcSet)
		strcpy(SOURCEPATH, "./");
	else if (SOURCEPATH[strlen(SOURCEPATH) - 1] != '/')
		strcat(SOURCEPATH, "/");
	printf("Source path: %s\n", SOURCEPATH);
	if (!dstSet)
		strcpy(DESTPATH, "./");
	else if (DESTPATH[strlen(DESTPATH) - 1] != '/')
		strcat(DESTPATH, "/");
	printf("Destination path: %s\n", DESTPATH);
	printf("Utilized cpu cores: %d\n", numCores);
	printf("Executing: ");
	if (phaseSet)
		printf("Phase #%d\n", phaseNum);
	else
		printf("All phases\n");
	printf("Processing from ");
	printf(FNAME, SLICE_START);
	printf(" (slice #:%d) to ", SLICE_START);
	printf(FNAME, SLICE_END);
	printf(" (slice #:%d)...\n", SLICE_END);

	if (phaseSet) {
		switch (phaseNum) {
		case 1:
			startPhase1();
			break;
		case 2:
			startPhase2();
			break;
		case 3:
			startPhase3();
			break;
		default:
			break;
		}
	} else {
		startPhase1();
		startPhase2();
		startPhase3();
	}
	printf(">>>> Completed!\n\n");

	return 0;
}
