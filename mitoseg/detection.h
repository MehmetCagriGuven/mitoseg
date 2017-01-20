/*
 * detection.h
 *
 *  Created on: Jan 15, 2017
 *      Author: fst
 */

#ifndef DETECTION_H_
#define DETECTION_H_

#include "opencv/cv.h"
#include "opencv/highgui.h"

#include "snake25d.h"

#define elem_type double
elem_type quick_select(elem_type arr[], int n);
int isValidMitoSlice_25d(snake25d *s, int slice, dataPacket *dp);

#endif /* DETECTION_H_ */
