/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * nn_batch_data.c
 *
 * Code generation for function 'nn_batch_data'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "nn_batch.h"
#include "nn_batch_data.h"

/* Variable Definitions */
emlrtCTX emlrtRootTLSGlobal = NULL;
const volatile char_T *emlrtBreakCheckR2012bFlagVar = NULL;
emlrtContext emlrtContextGlobal = { true,/* bFirstTime */
  false,                               /* bInitialized */
  131450U,                             /* fVersionInfo */
  NULL,                                /* fErrorFunction */
  "nn_batch",                          /* fFunctionName */
  NULL,                                /* fRTCallStack */
  false,                               /* bDebugMode */
  { 2045744189U, 2170104910U, 2743257031U, 4284093946U },/* fSigWrd */
  NULL                                 /* fSigMem */
};

emlrtRSInfo k_emlrtRSI = { 21,         /* lineNo */
  "eml_int_forloop_overflow_check",    /* fcnName */
  "C:\\MATLAB\\R2017a\\toolbox\\eml\\lib\\matlab\\eml\\eml_int_forloop_overflow_check.m"/* pathName */
};

emlrtRSInfo o_emlrtRSI = { 86,         /* lineNo */
  "dot",                               /* fcnName */
  "C:\\MATLAB\\R2017a\\toolbox\\eml\\lib\\matlab\\specfun\\dot.m"/* pathName */
};

emlrtRSInfo p_emlrtRSI = { 32,         /* lineNo */
  "xdotc",                             /* fcnName */
  "C:\\MATLAB\\R2017a\\toolbox\\eml\\eml\\+coder\\+internal\\+blas\\xdotc.m"/* pathName */
};

emlrtRSInfo q_emlrtRSI = { 49,         /* lineNo */
  "xdot",                              /* fcnName */
  "C:\\MATLAB\\R2017a\\toolbox\\eml\\eml\\+coder\\+internal\\+blas\\xdot.m"/* pathName */
};

emlrtRSInfo r_emlrtRSI = { 50,         /* lineNo */
  "xdot",                              /* fcnName */
  "C:\\MATLAB\\R2017a\\toolbox\\eml\\eml\\+coder\\+internal\\+blas\\xdot.m"/* pathName */
};

emlrtRSInfo s_emlrtRSI = { 51,         /* lineNo */
  "xdot",                              /* fcnName */
  "C:\\MATLAB\\R2017a\\toolbox\\eml\\eml\\+coder\\+internal\\+blas\\xdot.m"/* pathName */
};

emlrtRSInfo t_emlrtRSI = { 10,         /* lineNo */
  "int",                               /* fcnName */
  "C:\\MATLAB\\R2017a\\toolbox\\eml\\eml\\+coder\\+internal\\+blas\\int.m"/* pathName */
};

emlrtRSInfo cb_emlrtRSI = { 97,        /* lineNo */
  "bsxfun",                            /* fcnName */
  "C:\\MATLAB\\R2017a\\toolbox\\eml\\lib\\matlab\\elmat\\bsxfun.m"/* pathName */
};

emlrtRSInfo db_emlrtRSI = { 91,        /* lineNo */
  "bsxfun",                            /* fcnName */
  "C:\\MATLAB\\R2017a\\toolbox\\eml\\lib\\matlab\\elmat\\bsxfun.m"/* pathName */
};

emlrtRSInfo eb_emlrtRSI = { 87,        /* lineNo */
  "bsxfun",                            /* fcnName */
  "C:\\MATLAB\\R2017a\\toolbox\\eml\\lib\\matlab\\elmat\\bsxfun.m"/* pathName */
};

emlrtRSInfo fb_emlrtRSI = { 85,        /* lineNo */
  "bsxfun",                            /* fcnName */
  "C:\\MATLAB\\R2017a\\toolbox\\eml\\lib\\matlab\\elmat\\bsxfun.m"/* pathName */
};

emlrtRSInfo ic_emlrtRSI = { 443,       /* lineNo */
  "sortIdx",                           /* fcnName */
  "C:\\MATLAB\\R2017a\\toolbox\\eml\\eml\\+coder\\+internal\\sortIdx.m"/* pathName */
};

emlrtRTEInfo d_emlrtRTEI = { 75,       /* lineNo */
  1,                                   /* colNo */
  "bsxfun",                            /* fName */
  "C:\\MATLAB\\R2017a\\toolbox\\eml\\lib\\matlab\\elmat\\bsxfun.m"/* pName */
};

emlrtRTEInfo e_emlrtRTEI = { 76,       /* lineNo */
  1,                                   /* colNo */
  "bsxfun",                            /* fName */
  "C:\\MATLAB\\R2017a\\toolbox\\eml\\lib\\matlab\\elmat\\bsxfun.m"/* pName */
};

emlrtRTEInfo f_emlrtRTEI = { 95,       /* lineNo */
  5,                                   /* colNo */
  "bsxfun",                            /* fName */
  "C:\\MATLAB\\R2017a\\toolbox\\eml\\lib\\matlab\\elmat\\bsxfun.m"/* pName */
};

emlrtRTEInfo v_emlrtRTEI = { 21,       /* lineNo */
  15,                                  /* colNo */
  "bsxfun",                            /* fName */
  "C:\\MATLAB\\R2017a\\toolbox\\eml\\lib\\matlab\\elmat\\bsxfun.m"/* pName */
};

emlrtRTEInfo w_emlrtRTEI = { 19,       /* lineNo */
  15,                                  /* colNo */
  "bsxfun",                            /* fName */
  "C:\\MATLAB\\R2017a\\toolbox\\eml\\lib\\matlab\\elmat\\bsxfun.m"/* pName */
};

/* End of code generation (nn_batch_data.c) */
