/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * nn_batch.h
 *
 * Code generation for function 'nn_batch'
 *
 */

#ifndef NN_BATCH_H
#define NN_BATCH_H

/* Include files */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mwmathutil.h"
#include "tmwtypes.h"
#include "mex.h"
#include "emlrt.h"
#include "rtwtypes.h"
#include "nn_batch_types.h"

/* Function Declarations */
extern void nn_batch(const emlrtStack *sp, const emxArray_real_T *X, const
                     emxArray_real_T *Y, real_T k, emxArray_real_T *neighbors,
                     emxArray_real_T *distances);

#endif

/* End of code generation (nn_batch.h) */
