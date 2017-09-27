/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * error.c
 *
 * Code generation for function 'error'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "nn_batch.h"
#include "error.h"

/* Variable Definitions */
static emlrtMCInfo emlrtMCI = { 27,    /* lineNo */
  5,                                   /* colNo */
  "error",                             /* fName */
  "C:\\MATLAB\\R2017a\\toolbox\\eml\\lib\\matlab\\lang\\error.m"/* pName */
};

static emlrtRSInfo nc_emlrtRSI = { 27, /* lineNo */
  "error",                             /* fcnName */
  "C:\\MATLAB\\R2017a\\toolbox\\eml\\lib\\matlab\\lang\\error.m"/* pathName */
};

/* Function Declarations */
static void b_error(const emlrtStack *sp, const mxArray *b, emlrtMCInfo
                    *location);

/* Function Definitions */
static void b_error(const emlrtStack *sp, const mxArray *b, emlrtMCInfo
                    *location)
{
  const mxArray *pArray;
  pArray = b;
  emlrtCallMATLABR2012b(sp, 0, NULL, 1, &pArray, "error", true, location);
}

void error(const emlrtStack *sp)
{
  const mxArray *y;
  const mxArray *m0;
  static const int32_T iv1[2] = { 1, 46 };

  static const char_T varargin_1[46] = { 'r', 'e', 'q', 'u', 'e', 's', 't', 'e',
    'd', ' ', 'm', 'o', 'r', 'e', ' ', 'n', 'e', 'i', 'g', 'h', 'b', 'o', 'r',
    's', ' ', 't', 'h', 'a', 'n', ' ', 'n', 'u', 'm', 'b', 'e', 'r', ' ', 'o',
    'f', ' ', 'p', 'o', 'i', 'n', 't', 's' };

  emlrtStack st;
  st.prev = sp;
  st.tls = sp->tls;
  y = NULL;
  m0 = emlrtCreateCharArray(2, iv1);
  emlrtInitCharArrayR2013a(sp, 46, m0, &varargin_1[0]);
  emlrtAssign(&y, m0);
  st.site = &nc_emlrtRSI;
  b_error(&st, y, &emlrtMCI);
}

/* End of code generation (error.c) */
