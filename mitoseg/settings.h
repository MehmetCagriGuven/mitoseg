/*
 * settings.h
 *
 *  Created on: Jan 13, 2017
 *      Author: fst
 */

#ifndef SETTINGS_H_
#define SETTINGS_H_

///////////// PHASE 1 PARAMS /////////////
// Visualization params
#define VISUAL_LE_SSIZE 8

// Preprocessing params
#define TFACTOR 2.0							//nm/px
#define AUTOLEVELSCUT 0.005
#define BSMOOTH_SPAT (60 / TFACTOR)			//nm -> px
#define BSMOOTH_GRAY 0.2f
#define SMOOTH_STDEV (3.0 / TFACTOR)		//nm -> px

// Common curve fit params
#define FC_MAXCURV 1.0
#define FC_MAXITER 200
#define FC_STEP 1.0
#define FC_XYRANGE 2 // fitcurve2-3
#define FC_XYSTEP 1 // fitcurve2-3
#define FC_HRANGE 2 // fitcurve2-3
#define FC_HSTEP 1 // fitcurve2-3
#define FC_TOL 0.01 // fitcurve1
#define LE_BINS 4

// Low freq. curve fit params
#define LE_SSIZE_LO (4 / TFACTOR)								//nm -> px/blk
#define LE_BSIZE_LO (30 / TFACTOR)								//nm -> px/blk
#define FCS_XYSTEP_LO ((int)((60 / TFACTOR) / LE_SSIZE_LO))		//nm -> blk
#define FCS_XYRANGE_LO ((int)((60 / TFACTOR) / LE_SSIZE_LO))	//nm -> blk
#define FCS_MINLEN_LO ((int)((100 / TFACTOR) / LE_SSIZE_LO))	//nm -> blk
#define FLC_THRESH_LO (100.0 * LE_BSIZE_LO / LE_SSIZE_LO)		//pt -> pt (scaled)
#define FLC_AVGTHRESH_LO (1.2 * LE_BSIZE_LO)					//pt/px -> pt/blk
#define FCS_MAXCURVE_LO 20000
#define FCS_INIT_THRESH_LO 0.4f

// High freq. curve fit params
#define LE_SSIZE_HI (4 / TFACTOR)								//nm -> px/blk
#define LE_BSIZE_HI (8 / TFACTOR)								//nm -> px/blk
#define FCS_XYSTEP_HI ((int)((16 / TFACTOR) / LE_SSIZE_HI))		//nm -> blk
#define FCS_XYRANGE_HI ((int)((16 / TFACTOR) / LE_SSIZE_HI))	//nm -> blk
#define FCS_MINLEN_HI ((int)((20 / TFACTOR) / LE_SSIZE_HI))		//nm -> blk
#define FLC_THRESH_HI (20 * LE_BSIZE_HI / LE_SSIZE_HI)			//pt -> pt (scaled)
#define FLC_AVGTHRESH_HI (0.7 * LE_BSIZE_HI)					//pt/px -> pt/blk
#define FCS_MAXCURVE_HI 20000
#define FCS_INIT_THRESH_HI 0.4f

// Low Freq. curves coverage distance
#define CURV_HILO_COV ((int)((60 / TFACTOR) / LE_SSIZE_LO))		//nm -> blk
#define CURV_HILO_COV_TH 0.7

// Detection params
#define DT_SET				3

#define DT_BOUNDARY_T 20.0
#define DT_REGION_T 0.1

#if DT_SET == 1
// old settings for window 0 ayar 2
#define DT_GAP_T 5.0
#define DT_GAPTOTAL_T ((500.0 / TFACTOR) / LE_SSIZE_LO)		//nm -> blk
#define DT_GAPMAX_T ((400.0 / TFACTOR) / LE_SSIZE_LO)		//nm -> blk
#define DT_GAPTOTALRATIO_T 0.4
#define DT_GAPBORDERRATIO_T 0.3
#define DT_GAPBORDERCOUNT_T 1
#elif DT_SET == 2
// window 3 ayar 0
#define DT_GAP_T 10.0
#define DT_GAPTOTAL_T ((500.0 / TFACTOR) / LE_SSIZE_LO)		//nm -> blk
#define DT_GAPMAX_T ((400.0 / TFACTOR) / LE_SSIZE_LO)		//nm -> blk
#define DT_GAPTOTALRATIO_T 0.4
#define DT_GAPBORDERRATIO_T 0.3
#define DT_GAPBORDERCOUNT_T 1
#elif DT_SET == 3
// window 5 ayar 1
#define DT_GAP_T 15.0
#define DT_GAPTOTAL_T ((600.0 / TFACTOR) / LE_SSIZE_LO)		//nm -> blk
#define DT_GAPMAX_T ((600.0 / TFACTOR) / LE_SSIZE_LO)		//nm -> blk
#define DT_GAPTOTALRATIO_T 0.4
#define DT_GAPBORDERRATIO_T 0.4
#define DT_GAPBORDERCOUNT_T 1
#endif

#define DT_CURVG_T 15.0
#define DT_CURVL_T 0.6
#define DT_CURVL_SEG 0.25
#define DT_SIGN_SMOOTH 0.05
#define DT_SIGN_RATIO 4.0
#define DT_SIGN_MAXNUMCRIT_T 4
#define DT_MAJORAXIS_T ((2000.0 / TFACTOR) / LE_SSIZE_LO)	//nm -> blk
#define DT_MINORAXIS_T ((140.0 / TFACTOR) / LE_SSIZE_LO)	//nm -> blk
#define DT_MINAXIS_T ((70.0 / TFACTOR) / LE_SSIZE_LO)		//nm -> blk
#define DT_MINMINOR_RATIO_T 0.2 //0.35
#define DT_MINAREA SN_MINAREA
#define DT_MAXAREA SN_MAXAREA
#define DT_REPORT 0
#define DT25D_MED_RANGE		5	//0 //3 //5 //7 //9
//////////////////////////////////////////

///////////// PHASE 2 PARAMS /////////////
// Common snake params
#define SN_N 100
#define SN_GAUSSIAN 1.0	//blk
#define SN_MAXAREA ((700000 / (TFACTOR*TFACTOR)) / (LE_SSIZE_LO*LE_SSIZE_LO))	//nm^2 -> blk^2
#define SN_MINAREA ((20000 / (TFACTOR*TFACTOR)) / (LE_SSIZE_LO*LE_SSIZE_LO))	//nm^2 -> blk^2
// 2.5D Snake params
#define SN25D_W_TENSION 1.0
#define SN25D_W_CURVATURE 200.0
#define SN25D_W_ZTENSION 5.0
#define SN25D_W_ZCURVATURE 5.0
#define SN25D_W_ECURVE 1.0
#define SN25D_W_EINF_MIN 0.5
#define SN25D_W_EINF_MAX 3.0
#define SN25D_W_EINF_STEP 0.5
#define SN25D_T				500
extern int sn25d_t;
#define SN25D_K 1.0
#define SN25D_INITR 10.0
#define SN25D_MAXITER 400000
#define SN25D_SHORTTERM_CONV_ITER 10000
#define SN25D_SHORTTERM_CONV 2.0
#define SN25D_LONGTERM_CONV_ITER 40000
#define SN25D_LONGTERM_CONV 10.0
#define SN25D_MAXAREA ((700000 / (TFACTOR*TFACTOR)) / (LE_SSIZE_LO*LE_SSIZE_LO))	//nm^2 -> blk^2
#define SN25D_MINAREA ((1000 / (TFACTOR*TFACTOR)) / (LE_SSIZE_LO*LE_SSIZE_LO))	    //nm^2 -> blk^2 (Vanishing area)
#define SN25D_INF_CONV 0.95
#define SN25D_INITPTS_EPS ((100 / TFACTOR) / LE_SSIZE_LO) //((150 / TFACTOR) / LE_SSIZE_LO)	//nm -> blk
#define SN25D_INITPTS_MIN ((int)(1.5 * sn25d_t))
//////////////////////////////////////////

///////////// PHASE 3 PARAMS /////////////
extern double POLY_VALIDITY;	// default: 0.75
#define POLY_MERGE 0.3
//////////////////////////////////////////

/////////// VISUALIZER PARAMS ////////////
#define MAXIMAGES 150
#define MAXBUFFER 150
#define WINDOW_WIDTH 720
#define WINDOW_HEIGHT 720
//////////////////////////////////////////

#endif /* SETTINGS_H_ */
