/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * nn_batch_initialize.c
 *
 * Code generation for function 'nn_batch_initialize'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "nn_batch.h"
#include "nn_batch_initialize.h"
#include "_coder_nn_batch_mex.h"
#include "nn_batch_data.h"

/* Function Definitions */
void nn_batch_initialize(void)
{
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  mexFunctionCreateRootTLS();
  emlrtBreakCheckR2012bFlagVar = emlrtGetBreakCheckFlagAddressR2012b();
  st.tls = emlrtRootTLSGlobal;
  emlrtClearAllocCountR2012b(&st, false, 0U, 0);
  emlrtEnterRtStackR2012b(&st);
  emlrtFirstTimeR2012b(emlrtRootTLSGlobal);
}

/* End of code generation (nn_batch_initialize.c) */
