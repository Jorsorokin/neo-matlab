/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * find_blockIndex_range.h
 *
 * Code generation for function 'find_blockIndex_range'
 *
 */

#ifndef FIND_BLOCKINDEX_RANGE_H
#define FIND_BLOCKINDEX_RANGE_H

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
extern void find_blockIndex_range(const emlrtStack *sp, real_T n, real_T m,
  emxArray_real_T *start, emxArray_real_T *stop);

#endif

/* End of code generation (find_blockIndex_range.h) */
