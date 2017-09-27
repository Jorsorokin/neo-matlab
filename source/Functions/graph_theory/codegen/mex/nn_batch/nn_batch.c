/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * nn_batch.c
 *
 * Code generation for function 'nn_batch'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "nn_batch.h"
#include "nn_batch_emxutil.h"
#include "eml_int_forloop_overflow_check.h"
#include "sort1.h"
#include "bsxfun.h"
#include "dot.h"
#include "find_blockIndex_range.h"
#include "error.h"
#include "any.h"
#include "nn_batch_data.h"
#include "blas.h"

/* Variable Definitions */
static emlrtRSInfo emlrtRSI = { 6,     /* lineNo */
  "nn_batch",                          /* fcnName */
  "Z:\\matlabscripts\\SpikeSorting\\Functions\\graph_theory\\nn_batch.m"/* pathName */
};

static emlrtRSInfo b_emlrtRSI = { 12,  /* lineNo */
  "nn_batch",                          /* fcnName */
  "Z:\\matlabscripts\\SpikeSorting\\Functions\\graph_theory\\nn_batch.m"/* pathName */
};

static emlrtRSInfo c_emlrtRSI = { 18,  /* lineNo */
  "nn_batch",                          /* fcnName */
  "Z:\\matlabscripts\\SpikeSorting\\Functions\\graph_theory\\nn_batch.m"/* pathName */
};

static emlrtRSInfo d_emlrtRSI = { 19,  /* lineNo */
  "nn_batch",                          /* fcnName */
  "Z:\\matlabscripts\\SpikeSorting\\Functions\\graph_theory\\nn_batch.m"/* pathName */
};

static emlrtRSInfo e_emlrtRSI = { 22,  /* lineNo */
  "nn_batch",                          /* fcnName */
  "Z:\\matlabscripts\\SpikeSorting\\Functions\\graph_theory\\nn_batch.m"/* pathName */
};

static emlrtRSInfo u_emlrtRSI = { 68,  /* lineNo */
  "eml_mtimes_helper",                 /* fcnName */
  "C:\\MATLAB\\R2017a\\toolbox\\eml\\lib\\matlab\\ops\\eml_mtimes_helper.m"/* pathName */
};

static emlrtRSInfo v_emlrtRSI = { 21,  /* lineNo */
  "eml_mtimes_helper",                 /* fcnName */
  "C:\\MATLAB\\R2017a\\toolbox\\eml\\lib\\matlab\\ops\\eml_mtimes_helper.m"/* pathName */
};

static emlrtRSInfo w_emlrtRSI = { 87,  /* lineNo */
  "xgemm",                             /* fcnName */
  "C:\\MATLAB\\R2017a\\toolbox\\eml\\eml\\+coder\\+internal\\+blas\\xgemm.m"/* pathName */
};

static emlrtRSInfo x_emlrtRSI = { 89,  /* lineNo */
  "xgemm",                             /* fcnName */
  "C:\\MATLAB\\R2017a\\toolbox\\eml\\eml\\+coder\\+internal\\+blas\\xgemm.m"/* pathName */
};

static emlrtRSInfo y_emlrtRSI = { 90,  /* lineNo */
  "xgemm",                             /* fcnName */
  "C:\\MATLAB\\R2017a\\toolbox\\eml\\eml\\+coder\\+internal\\+blas\\xgemm.m"/* pathName */
};

static emlrtRSInfo ab_emlrtRSI = { 92, /* lineNo */
  "xgemm",                             /* fcnName */
  "C:\\MATLAB\\R2017a\\toolbox\\eml\\eml\\+coder\\+internal\\+blas\\xgemm.m"/* pathName */
};

static emlrtRSInfo bb_emlrtRSI = { 85, /* lineNo */
  "xgemm",                             /* fcnName */
  "C:\\MATLAB\\R2017a\\toolbox\\eml\\eml\\+coder\\+internal\\+blas\\xgemm.m"/* pathName */
};

static emlrtRSInfo gb_emlrtRSI = { 16, /* lineNo */
  "abs",                               /* fcnName */
  "C:\\MATLAB\\R2017a\\toolbox\\eml\\lib\\matlab\\elfun\\abs.m"/* pathName */
};

static emlrtRSInfo hb_emlrtRSI = { 67, /* lineNo */
  "applyScalarFunction",               /* fcnName */
  "C:\\MATLAB\\R2017a\\toolbox\\eml\\eml\\+coder\\+internal\\applyScalarFunction.m"/* pathName */
};

static emlrtRSInfo ib_emlrtRSI = { 15, /* lineNo */
  "sqrt",                              /* fcnName */
  "C:\\MATLAB\\R2017a\\toolbox\\eml\\lib\\matlab\\elfun\\sqrt.m"/* pathName */
};

static emlrtRSInfo jb_emlrtRSI = { 26, /* lineNo */
  "applyScalarFunctionInPlace",        /* fcnName */
  "C:\\MATLAB\\R2017a\\toolbox\\eml\\eml\\+coder\\+internal\\applyScalarFunctionInPlace.m"/* pathName */
};

static emlrtRSInfo kb_emlrtRSI = { 26, /* lineNo */
  "sort",                              /* fcnName */
  "C:\\MATLAB\\R2017a\\toolbox\\eml\\lib\\matlab\\datafun\\sort.m"/* pathName */
};

static emlrtRTEInfo emlrtRTEI = { 1,   /* lineNo */
  34,                                  /* colNo */
  "nn_batch",                          /* fName */
  "Z:\\matlabscripts\\SpikeSorting\\Functions\\graph_theory\\nn_batch.m"/* pName */
};

static emlrtRTEInfo b_emlrtRTEI = { 17,/* lineNo */
  5,                                   /* colNo */
  "nn_batch",                          /* fName */
  "Z:\\matlabscripts\\SpikeSorting\\Functions\\graph_theory\\nn_batch.m"/* pName */
};

static emlrtRTEInfo c_emlrtRTEI = { 18,/* lineNo */
  5,                                   /* colNo */
  "nn_batch",                          /* fName */
  "Z:\\matlabscripts\\SpikeSorting\\Functions\\graph_theory\\nn_batch.m"/* pName */
};

static emlrtBCInfo emlrtBCI = { -1,    /* iFirst */
  -1,                                  /* iLast */
  17,                                  /* lineNo */
  19,                                  /* colNo */
  "start",                             /* aName */
  "nn_batch",                          /* fName */
  "Z:\\matlabscripts\\SpikeSorting\\Functions\\graph_theory\\nn_batch.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo b_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  17,                                  /* lineNo */
  27,                                  /* colNo */
  "stop",                              /* aName */
  "nn_batch",                          /* fName */
  "Z:\\matlabscripts\\SpikeSorting\\Functions\\graph_theory\\nn_batch.m",/* pName */
  0                                    /* checkKind */
};

static emlrtDCInfo emlrtDCI = { 17,    /* lineNo */
  13,                                  /* colNo */
  "nn_batch",                          /* fName */
  "Z:\\matlabscripts\\SpikeSorting\\Functions\\graph_theory\\nn_batch.m",/* pName */
  1                                    /* checkKind */
};

static emlrtBCInfo c_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  17,                                  /* lineNo */
  13,                                  /* colNo */
  "X",                                 /* aName */
  "nn_batch",                          /* fName */
  "Z:\\matlabscripts\\SpikeSorting\\Functions\\graph_theory\\nn_batch.m",/* pName */
  0                                    /* checkKind */
};

static emlrtDCInfo b_emlrtDCI = { 17,  /* lineNo */
  22,                                  /* colNo */
  "nn_batch",                          /* fName */
  "Z:\\matlabscripts\\SpikeSorting\\Functions\\graph_theory\\nn_batch.m",/* pName */
  1                                    /* checkKind */
};

static emlrtBCInfo d_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  17,                                  /* lineNo */
  22,                                  /* colNo */
  "X",                                 /* aName */
  "nn_batch",                          /* fName */
  "Z:\\matlabscripts\\SpikeSorting\\Functions\\graph_theory\\nn_batch.m",/* pName */
  0                                    /* checkKind */
};

static emlrtECInfo emlrtECI = { -1,    /* nDims */
  22,                                  /* lineNo */
  6,                                   /* colNo */
  "nn_batch",                          /* fName */
  "Z:\\matlabscripts\\SpikeSorting\\Functions\\graph_theory\\nn_batch.m"/* pName */
};

static emlrtBCInfo e_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  23,                                  /* lineNo */
  41,                                  /* colNo */
  "idx",                               /* aName */
  "nn_batch",                          /* fName */
  "Z:\\matlabscripts\\SpikeSorting\\Functions\\graph_theory\\nn_batch.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo f_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  23,                                  /* lineNo */
  43,                                  /* colNo */
  "idx",                               /* aName */
  "nn_batch",                          /* fName */
  "Z:\\matlabscripts\\SpikeSorting\\Functions\\graph_theory\\nn_batch.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo g_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  24,                                  /* lineNo */
  39,                                  /* colNo */
  "S",                                 /* aName */
  "nn_batch",                          /* fName */
  "Z:\\matlabscripts\\SpikeSorting\\Functions\\graph_theory\\nn_batch.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo h_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  24,                                  /* lineNo */
  41,                                  /* colNo */
  "S",                                 /* aName */
  "nn_batch",                          /* fName */
  "Z:\\matlabscripts\\SpikeSorting\\Functions\\graph_theory\\nn_batch.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo i_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  23,                                  /* lineNo */
  21,                                  /* colNo */
  "start",                             /* aName */
  "nn_batch",                          /* fName */
  "Z:\\matlabscripts\\SpikeSorting\\Functions\\graph_theory\\nn_batch.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo j_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  23,                                  /* lineNo */
  29,                                  /* colNo */
  "stop",                              /* aName */
  "nn_batch",                          /* fName */
  "Z:\\matlabscripts\\SpikeSorting\\Functions\\graph_theory\\nn_batch.m",/* pName */
  0                                    /* checkKind */
};

static emlrtDCInfo c_emlrtDCI = { 23,  /* lineNo */
  15,                                  /* colNo */
  "nn_batch",                          /* fName */
  "Z:\\matlabscripts\\SpikeSorting\\Functions\\graph_theory\\nn_batch.m",/* pName */
  1                                    /* checkKind */
};

static emlrtBCInfo k_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  23,                                  /* lineNo */
  15,                                  /* colNo */
  "neighbors",                         /* aName */
  "nn_batch",                          /* fName */
  "Z:\\matlabscripts\\SpikeSorting\\Functions\\graph_theory\\nn_batch.m",/* pName */
  0                                    /* checkKind */
};

static emlrtDCInfo d_emlrtDCI = { 23,  /* lineNo */
  24,                                  /* colNo */
  "nn_batch",                          /* fName */
  "Z:\\matlabscripts\\SpikeSorting\\Functions\\graph_theory\\nn_batch.m",/* pName */
  1                                    /* checkKind */
};

static emlrtBCInfo l_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  23,                                  /* lineNo */
  24,                                  /* colNo */
  "neighbors",                         /* aName */
  "nn_batch",                          /* fName */
  "Z:\\matlabscripts\\SpikeSorting\\Functions\\graph_theory\\nn_batch.m",/* pName */
  0                                    /* checkKind */
};

static emlrtECInfo b_emlrtECI = { -1,  /* nDims */
  23,                                  /* lineNo */
  5,                                   /* colNo */
  "nn_batch",                          /* fName */
  "Z:\\matlabscripts\\SpikeSorting\\Functions\\graph_theory\\nn_batch.m"/* pName */
};

static emlrtBCInfo m_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  24,                                  /* lineNo */
  21,                                  /* colNo */
  "start",                             /* aName */
  "nn_batch",                          /* fName */
  "Z:\\matlabscripts\\SpikeSorting\\Functions\\graph_theory\\nn_batch.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo n_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  24,                                  /* lineNo */
  29,                                  /* colNo */
  "stop",                              /* aName */
  "nn_batch",                          /* fName */
  "Z:\\matlabscripts\\SpikeSorting\\Functions\\graph_theory\\nn_batch.m",/* pName */
  0                                    /* checkKind */
};

static emlrtDCInfo e_emlrtDCI = { 24,  /* lineNo */
  15,                                  /* colNo */
  "nn_batch",                          /* fName */
  "Z:\\matlabscripts\\SpikeSorting\\Functions\\graph_theory\\nn_batch.m",/* pName */
  1                                    /* checkKind */
};

static emlrtBCInfo o_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  24,                                  /* lineNo */
  15,                                  /* colNo */
  "distances",                         /* aName */
  "nn_batch",                          /* fName */
  "Z:\\matlabscripts\\SpikeSorting\\Functions\\graph_theory\\nn_batch.m",/* pName */
  0                                    /* checkKind */
};

static emlrtDCInfo f_emlrtDCI = { 24,  /* lineNo */
  24,                                  /* colNo */
  "nn_batch",                          /* fName */
  "Z:\\matlabscripts\\SpikeSorting\\Functions\\graph_theory\\nn_batch.m",/* pName */
  1                                    /* checkKind */
};

static emlrtBCInfo p_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  24,                                  /* lineNo */
  24,                                  /* colNo */
  "distances",                         /* aName */
  "nn_batch",                          /* fName */
  "Z:\\matlabscripts\\SpikeSorting\\Functions\\graph_theory\\nn_batch.m",/* pName */
  0                                    /* checkKind */
};

static emlrtECInfo c_emlrtECI = { -1,  /* nDims */
  24,                                  /* lineNo */
  5,                                   /* colNo */
  "nn_batch",                          /* fName */
  "Z:\\matlabscripts\\SpikeSorting\\Functions\\graph_theory\\nn_batch.m"/* pName */
};

static emlrtRTEInfo t_emlrtRTEI = { 104,/* lineNo */
  23,                                  /* colNo */
  "eml_mtimes_helper",                 /* fName */
  "C:\\MATLAB\\R2017a\\toolbox\\eml\\lib\\matlab\\ops\\eml_mtimes_helper.m"/* pName */
};

static emlrtRTEInfo u_emlrtRTEI = { 99,/* lineNo */
  23,                                  /* colNo */
  "eml_mtimes_helper",                 /* fName */
  "C:\\MATLAB\\R2017a\\toolbox\\eml\\lib\\matlab\\ops\\eml_mtimes_helper.m"/* pName */
};

static emlrtDCInfo g_emlrtDCI = { 10,  /* lineNo */
  22,                                  /* colNo */
  "nn_batch",                          /* fName */
  "Z:\\matlabscripts\\SpikeSorting\\Functions\\graph_theory\\nn_batch.m",/* pName */
  1                                    /* checkKind */
};

static emlrtDCInfo h_emlrtDCI = { 10,  /* lineNo */
  22,                                  /* colNo */
  "nn_batch",                          /* fName */
  "Z:\\matlabscripts\\SpikeSorting\\Functions\\graph_theory\\nn_batch.m",/* pName */
  4                                    /* checkKind */
};

static emlrtDCInfo i_emlrtDCI = { 11,  /* lineNo */
  21,                                  /* colNo */
  "nn_batch",                          /* fName */
  "Z:\\matlabscripts\\SpikeSorting\\Functions\\graph_theory\\nn_batch.m",/* pName */
  1                                    /* checkKind */
};

/* Function Definitions */
void nn_batch(const emlrtStack *sp, const emxArray_real_T *X, const
              emxArray_real_T *Y, real_T k, emxArray_real_T *neighbors,
              emxArray_real_T *distances)
{
  boolean_T b_k[2];
  int32_T i0;
  int32_T loop_ub;
  emxArray_real_T *start;
  emxArray_real_T *stop;
  int32_T i;
  emxArray_real_T *pts;
  emxArray_real_T *S;
  emxArray_int32_T *idx;
  emxArray_int32_T *r0;
  emxArray_int32_T *r1;
  emxArray_real_T *b;
  emxArray_real_T *a;
  emxArray_real_T *av;
  emxArray_real_T *bv;
  emxArray_real_T *cv;
  emxArray_int32_T *iidx;
  emxArray_real_T *r2;
  int32_T b_b;
  real_T alpha1;
  int32_T ck;
  int32_T csz_idx_1;
  int32_T nc1;
  uint32_T uv0[2];
  int32_T bsub;
  int32_T na1;
  int32_T nb1;
  boolean_T overflow;
  uint32_T a_idx_0;
  real_T beta1;
  char_T TRANSB;
  char_T TRANSA;
  ptrdiff_t m_t;
  ptrdiff_t n_t;
  ptrdiff_t k_t;
  ptrdiff_t lda_t;
  ptrdiff_t ldb_t;
  ptrdiff_t ldc_t;
  int32_T iv0[2];
  emlrtStack st;
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack d_st;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  d_st.prev = &c_st;
  d_st.tls = c_st.tls;
  emlrtHeapReferenceStackEnterFcnR2012b(sp);
  b_k[0] = (k > (uint32_T)X->size[0]);
  b_k[1] = (k > (uint32_T)Y->size[0]);
  if (any(b_k)) {
    st.site = &emlrtRSI;
    error(&st);
  }

  /*  calculate the ranges for looping */
  i0 = neighbors->size[0] * neighbors->size[1];
  neighbors->size[0] = X->size[0];
  if (!(k >= 0.0)) {
    emlrtNonNegativeCheckR2012b(k, &h_emlrtDCI, sp);
  }

  if (k != (int32_T)muDoubleScalarFloor(k)) {
    emlrtIntegerCheckR2012b(k, &g_emlrtDCI, sp);
  }

  neighbors->size[1] = (int32_T)k;
  emxEnsureCapacity(sp, (emxArray__common *)neighbors, i0, sizeof(real_T),
                    &emlrtRTEI);
  if (!(k >= 0.0)) {
    emlrtNonNegativeCheckR2012b(k, &h_emlrtDCI, sp);
  }

  if (k != (int32_T)muDoubleScalarFloor(k)) {
    emlrtIntegerCheckR2012b(k, &g_emlrtDCI, sp);
  }

  loop_ub = X->size[0] * (int32_T)k;
  for (i0 = 0; i0 < loop_ub; i0++) {
    neighbors->data[i0] = 0.0;
  }

  i0 = distances->size[0] * distances->size[1];
  distances->size[0] = X->size[0];
  if (k != (int32_T)muDoubleScalarFloor(k)) {
    emlrtIntegerCheckR2012b(k, &i_emlrtDCI, sp);
  }

  distances->size[1] = (int32_T)k;
  emxEnsureCapacity(sp, (emxArray__common *)distances, i0, sizeof(real_T),
                    &emlrtRTEI);
  if (k != (int32_T)muDoubleScalarFloor(k)) {
    emlrtIntegerCheckR2012b(k, &i_emlrtDCI, sp);
  }

  loop_ub = X->size[0] * (int32_T)k;
  for (i0 = 0; i0 < loop_ub; i0++) {
    distances->data[i0] = 0.0;
  }

  emxInit_real_T(sp, &start, 2, &emlrtRTEI, true);
  emxInit_real_T(sp, &stop, 2, &emlrtRTEI, true);
  st.site = &b_emlrtRSI;
  find_blockIndex_range(&st, Y->size[0], X->size[0], start, stop);
  i = 0;
  emxInit_real_T(sp, &pts, 2, &b_emlrtRTEI, true);
  emxInit_real_T(sp, &S, 2, &c_emlrtRTEI, true);
  emxInit_int32_T(sp, &idx, 2, &emlrtRTEI, true);
  emxInit_int32_T1(sp, &r0, 1, &emlrtRTEI, true);
  emxInit_int32_T1(sp, &r1, 1, &emlrtRTEI, true);
  emxInit_real_T(sp, &b, 2, &emlrtRTEI, true);
  emxInit_real_T1(sp, &a, 1, &emlrtRTEI, true);
  emxInit_real_T1(sp, &av, 1, &d_emlrtRTEI, true);
  emxInit_real_T1(sp, &bv, 1, &e_emlrtRTEI, true);
  emxInit_real_T1(sp, &cv, 1, &f_emlrtRTEI, true);
  emxInit_int32_T(sp, &iidx, 2, &emlrtRTEI, true);
  emxInit_real_T(sp, &r2, 2, &emlrtRTEI, true);
  while (i <= start->size[1] - 1) {
    /*  euclidean distance */
    i0 = start->size[1];
    b_b = i + 1;
    if (!((b_b >= 1) && (b_b <= i0))) {
      emlrtDynamicBoundsCheckR2012b(b_b, 1, i0, &emlrtBCI, sp);
    }

    i0 = stop->size[1];
    b_b = i + 1;
    if (!((b_b >= 1) && (b_b <= i0))) {
      emlrtDynamicBoundsCheckR2012b(b_b, 1, i0, &b_emlrtBCI, sp);
    }

    if (start->data[i] > stop->data[i]) {
      b_b = 1;
      i0 = 1;
    } else {
      alpha1 = start->data[i];
      if (alpha1 != (int32_T)muDoubleScalarFloor(alpha1)) {
        emlrtIntegerCheckR2012b(alpha1, &emlrtDCI, sp);
      }

      i0 = X->size[0];
      b_b = (int32_T)alpha1;
      if (!((b_b >= 1) && (b_b <= i0))) {
        emlrtDynamicBoundsCheckR2012b(b_b, 1, i0, &c_emlrtBCI, sp);
      }

      alpha1 = stop->data[i];
      if (alpha1 != (int32_T)muDoubleScalarFloor(alpha1)) {
        emlrtIntegerCheckR2012b(alpha1, &b_emlrtDCI, sp);
      }

      i0 = X->size[0];
      ck = (int32_T)alpha1;
      if (!((ck >= 1) && (ck <= i0))) {
        emlrtDynamicBoundsCheckR2012b(ck, 1, i0, &d_emlrtBCI, sp);
      }

      i0 = ck + 1;
    }

    loop_ub = X->size[1];
    ck = pts->size[0] * pts->size[1];
    pts->size[0] = i0 - b_b;
    pts->size[1] = loop_ub;
    emxEnsureCapacity(sp, (emxArray__common *)pts, ck, sizeof(real_T),
                      &emlrtRTEI);
    for (ck = 0; ck < loop_ub; ck++) {
      csz_idx_1 = i0 - b_b;
      for (nc1 = 0; nc1 < csz_idx_1; nc1++) {
        pts->data[nc1 + pts->size[0] * ck] = X->data[((b_b + nc1) + X->size[0] *
          ck) - 1];
      }
    }

    st.site = &d_emlrtRSI;
    dot(&st, pts, Y, av);
    st.site = &d_emlrtRSI;
    ck = b->size[0] * b->size[1];
    b->size[0] = Y->size[1];
    b->size[1] = Y->size[0];
    emxEnsureCapacity(&st, (emxArray__common *)b, ck, sizeof(real_T), &emlrtRTEI);
    loop_ub = Y->size[0];
    for (ck = 0; ck < loop_ub; ck++) {
      csz_idx_1 = Y->size[1];
      for (nc1 = 0; nc1 < csz_idx_1; nc1++) {
        b->data[nc1 + b->size[0] * ck] = Y->data[ck + Y->size[0] * nc1];
      }
    }

    b_st.site = &v_emlrtRSI;
    ck = X->size[1];
    if (!(ck == b->size[0])) {
      ck = X->size[1];
      if (((i0 - b_b == 1) && (ck == 1)) || ((b->size[0] == 1) && (b->size[1] ==
            1))) {
        emlrtErrorWithMessageIdR2012b(&b_st, &u_emlrtRTEI,
          "Coder:toolbox:mtimes_noDynamicScalarExpansion", 0);
      } else {
        emlrtErrorWithMessageIdR2012b(&b_st, &t_emlrtRTEI,
          "Coder:MATLAB:innerdim", 0);
      }
    }

    ck = X->size[1];
    if ((ck == 1) || (b->size[0] == 1)) {
      i0 = S->size[0] * S->size[1];
      S->size[0] = pts->size[0];
      S->size[1] = b->size[1];
      emxEnsureCapacity(&st, (emxArray__common *)S, i0, sizeof(real_T),
                        &emlrtRTEI);
      loop_ub = pts->size[0];
      for (i0 = 0; i0 < loop_ub; i0++) {
        csz_idx_1 = b->size[1];
        for (b_b = 0; b_b < csz_idx_1; b_b++) {
          S->data[i0 + S->size[0] * b_b] = 0.0;
          bsub = pts->size[1];
          for (ck = 0; ck < bsub; ck++) {
            S->data[i0 + S->size[0] * b_b] += pts->data[i0 + pts->size[0] * ck] *
              b->data[ck + b->size[0] * b_b];
          }
        }
      }
    } else {
      uv0[1] = (uint32_T)b->size[1];
      ck = S->size[0] * S->size[1];
      S->size[0] = i0 - b_b;
      S->size[1] = (int32_T)uv0[1];
      emxEnsureCapacity(&st, (emxArray__common *)S, ck, sizeof(real_T),
                        &emlrtRTEI);
      loop_ub = S->size[1];
      for (ck = 0; ck < loop_ub; ck++) {
        csz_idx_1 = S->size[0];
        for (nc1 = 0; nc1 < csz_idx_1; nc1++) {
          S->data[nc1 + S->size[0] * ck] = 0.0;
        }
      }

      b_st.site = &u_emlrtRSI;
      if ((i0 - b_b < 1) || (b->size[1] < 1)) {
      } else {
        ck = X->size[1];
        if (ck < 1) {
        } else {
          c_st.site = &w_emlrtRSI;
          c_st.site = &w_emlrtRSI;
          c_st.site = &w_emlrtRSI;
          c_st.site = &x_emlrtRSI;
          c_st.site = &y_emlrtRSI;
          c_st.site = &ab_emlrtRSI;
          c_st.site = &bb_emlrtRSI;
          alpha1 = 1.0;
          beta1 = 0.0;
          TRANSB = 'N';
          TRANSA = 'N';
          m_t = (ptrdiff_t)(i0 - b_b);
          n_t = (ptrdiff_t)b->size[1];
          ck = X->size[1];
          k_t = (ptrdiff_t)ck;
          lda_t = (ptrdiff_t)(i0 - b_b);
          ck = X->size[1];
          ldb_t = (ptrdiff_t)ck;
          ldc_t = (ptrdiff_t)(i0 - b_b);
          dgemm(&TRANSA, &TRANSB, &m_t, &n_t, &k_t, &alpha1, &pts->data[0],
                &lda_t, &b->data[0], &ldb_t, &beta1, &S->data[0], &ldc_t);
        }
      }
    }

    st.site = &c_emlrtRSI;
    dot(&st, pts, Y, a);
    i0 = r2->size[0] * r2->size[1];
    r2->size[0] = S->size[0];
    r2->size[1] = S->size[1];
    emxEnsureCapacity(sp, (emxArray__common *)r2, i0, sizeof(real_T), &emlrtRTEI);
    loop_ub = S->size[0] * S->size[1];
    for (i0 = 0; i0 < loop_ub; i0++) {
      r2->data[i0] = 2.0 * S->data[i0];
    }

    st.site = &d_emlrtRSI;
    bsxfun(&st, av, r2, b);
    st.site = &c_emlrtRSI;
    na1 = a->size[0];
    nb1 = b->size[0];
    if ((a->size[0] != 1) && (b->size[0] != 1) && (a->size[0] != b->size[0])) {
      overflow = false;
    } else {
      overflow = true;
    }

    if (!overflow) {
      emlrtErrorWithMessageIdR2012b(&st, &w_emlrtRTEI,
        "MATLAB:bsxfun:arrayDimensionsMustMatch", 0);
    }

    overflow = !(a->size[0] != b->size[0]);
    if (!overflow) {
      emlrtErrorWithMessageIdR2012b(&st, &v_emlrtRTEI,
        "Coder:toolbox:bsxfun_dynamicExpansion", 0);
    }

    if (a->size[0] <= b->size[0]) {
      bsub = a->size[0];
    } else {
      bsub = b->size[0];
    }

    csz_idx_1 = b->size[1];
    i0 = pts->size[0] * pts->size[1];
    pts->size[0] = bsub;
    pts->size[1] = csz_idx_1;
    emxEnsureCapacity(&st, (emxArray__common *)pts, i0, sizeof(real_T),
                      &emlrtRTEI);
    if (!((pts->size[0] == 0) || (pts->size[1] == 0))) {
      a_idx_0 = (uint32_T)a->size[0];
      i0 = av->size[0];
      av->size[0] = (int32_T)a_idx_0;
      emxEnsureCapacity(&st, (emxArray__common *)av, i0, sizeof(real_T),
                        &emlrtRTEI);
      a_idx_0 = (uint32_T)b->size[0];
      i0 = bv->size[0];
      bv->size[0] = (int32_T)a_idx_0;
      emxEnsureCapacity(&st, (emxArray__common *)bv, i0, sizeof(real_T),
                        &emlrtRTEI);
      bsub = 1;
      csz_idx_1 = 0;
      nc1 = pts->size[0];
      b_b = pts->size[0] * pts->size[1] - pts->size[0];
      b_st.site = &fb_emlrtRSI;
      overflow = ((!(0 > b_b)) && (b_b > MAX_int32_T - pts->size[0]));
      if (overflow) {
        c_st.site = &k_emlrtRSI;
        check_forloop_overflow_error(&c_st);
      }

      for (ck = 0; ck <= b_b; ck += nc1) {
        b_st.site = &eb_emlrtRSI;
        if ((!(1 > na1)) && (na1 > 2147483646)) {
          c_st.site = &k_emlrtRSI;
          check_forloop_overflow_error(&c_st);
        }

        for (loop_ub = 0; loop_ub + 1 <= na1; loop_ub++) {
          av->data[loop_ub] = a->data[loop_ub];
        }

        b_st.site = &db_emlrtRSI;
        if ((!(1 > nb1)) && (nb1 > 2147483646)) {
          c_st.site = &k_emlrtRSI;
          check_forloop_overflow_error(&c_st);
        }

        for (loop_ub = 0; loop_ub + 1 <= nb1; loop_ub++) {
          bv->data[loop_ub] = b->data[csz_idx_1 + loop_ub];
        }

        i0 = cv->size[0];
        cv->size[0] = av->size[0];
        emxEnsureCapacity(&st, (emxArray__common *)cv, i0, sizeof(real_T),
                          &emlrtRTEI);
        loop_ub = av->size[0];
        for (i0 = 0; i0 < loop_ub; i0++) {
          cv->data[i0] = av->data[i0] + bv->data[i0];
        }

        b_st.site = &cb_emlrtRSI;
        if (nc1 > 2147483646) {
          c_st.site = &k_emlrtRSI;
          check_forloop_overflow_error(&c_st);
        }

        for (loop_ub = 0; loop_ub + 1 <= nc1; loop_ub++) {
          pts->data[ck + loop_ub] = cv->data[loop_ub];
        }

        if (bsub < b->size[1]) {
          csz_idx_1 += nb1;
          bsub++;
        } else {
          bsub = 1;
        }
      }
    }

    st.site = &c_emlrtRSI;
    b_st.site = &gb_emlrtRSI;
    for (i0 = 0; i0 < 2; i0++) {
      uv0[i0] = (uint32_T)pts->size[i0];
    }

    i0 = S->size[0] * S->size[1];
    S->size[0] = (int32_T)uv0[0];
    S->size[1] = (int32_T)uv0[1];
    emxEnsureCapacity(&b_st, (emxArray__common *)S, i0, sizeof(real_T),
                      &emlrtRTEI);
    bsub = pts->size[0] * pts->size[1];
    c_st.site = &hb_emlrtRSI;
    if ((!(1 > bsub)) && (bsub > 2147483646)) {
      d_st.site = &k_emlrtRSI;
      check_forloop_overflow_error(&d_st);
    }

    for (loop_ub = 0; loop_ub + 1 <= bsub; loop_ub++) {
      S->data[loop_ub] = muDoubleScalarAbs(pts->data[loop_ub]);
    }

    st.site = &c_emlrtRSI;
    b_st.site = &ib_emlrtRSI;
    bsub = S->size[0] * S->size[1];
    c_st.site = &jb_emlrtRSI;
    if ((!(1 > bsub)) && (bsub > 2147483646)) {
      d_st.site = &k_emlrtRSI;
      check_forloop_overflow_error(&d_st);
    }

    for (loop_ub = 0; loop_ub + 1 <= bsub; loop_ub++) {
      S->data[loop_ub] = muDoubleScalarSqrt(S->data[loop_ub]);
    }

    /*  nearest neighbors */
    st.site = &e_emlrtRSI;
    i0 = b->size[0] * b->size[1];
    b->size[0] = S->size[0];
    b->size[1] = S->size[1];
    emxEnsureCapacity(&st, (emxArray__common *)b, i0, sizeof(real_T), &emlrtRTEI);
    loop_ub = S->size[0] * S->size[1];
    for (i0 = 0; i0 < loop_ub; i0++) {
      b->data[i0] = S->data[i0];
    }

    b_st.site = &kb_emlrtRSI;
    sort(&b_st, b, iidx);
    i0 = idx->size[0] * idx->size[1];
    idx->size[0] = iidx->size[0];
    idx->size[1] = iidx->size[1];
    emxEnsureCapacity(&st, (emxArray__common *)idx, i0, sizeof(int32_T),
                      &emlrtRTEI);
    loop_ub = iidx->size[0] * iidx->size[1];
    for (i0 = 0; i0 < loop_ub; i0++) {
      idx->data[i0] = iidx->data[i0];
    }

    i0 = S->size[0] * S->size[1];
    b_b = b->size[0] * b->size[1];
    if (i0 != b_b) {
      emlrtSizeEqCheck1DR2012b(i0, b_b, &emlrtECI, sp);
    }

    /*  sorts smallest - largest columns for each row */
    if (2.0 > k + 1.0) {
      i0 = 0;
      ck = 0;
    } else {
      i0 = idx->size[0];
      if (!(2 <= i0)) {
        emlrtDynamicBoundsCheckR2012b(2, 1, i0, &e_emlrtBCI, sp);
      }

      i0 = 1;
      b_b = idx->size[0];
      ck = (int32_T)(k + 1.0);
      if (!((ck >= 1) && (ck <= b_b))) {
        emlrtDynamicBoundsCheckR2012b(ck, 1, b_b, &f_emlrtBCI, sp);
      }
    }

    b_b = start->size[1];
    nc1 = i + 1;
    if (!((nc1 >= 1) && (nc1 <= b_b))) {
      emlrtDynamicBoundsCheckR2012b(nc1, 1, b_b, &i_emlrtBCI, sp);
    }

    b_b = stop->size[1];
    nc1 = i + 1;
    if (!((nc1 >= 1) && (nc1 <= b_b))) {
      emlrtDynamicBoundsCheckR2012b(nc1, 1, b_b, &j_emlrtBCI, sp);
    }

    if (start->data[i] > stop->data[i]) {
      b_b = 0;
      bsub = 0;
    } else {
      alpha1 = start->data[i];
      if (alpha1 != (int32_T)muDoubleScalarFloor(alpha1)) {
        emlrtIntegerCheckR2012b(alpha1, &c_emlrtDCI, sp);
      }

      b_b = neighbors->size[0];
      nc1 = (int32_T)alpha1;
      if (!((nc1 >= 1) && (nc1 <= b_b))) {
        emlrtDynamicBoundsCheckR2012b(nc1, 1, b_b, &k_emlrtBCI, sp);
      }

      b_b = nc1 - 1;
      alpha1 = stop->data[i];
      if (alpha1 != (int32_T)muDoubleScalarFloor(alpha1)) {
        emlrtIntegerCheckR2012b(alpha1, &d_emlrtDCI, sp);
      }

      nc1 = neighbors->size[0];
      bsub = (int32_T)alpha1;
      if (!((bsub >= 1) && (bsub <= nc1))) {
        emlrtDynamicBoundsCheckR2012b(bsub, 1, nc1, &l_emlrtBCI, sp);
      }
    }

    nc1 = r0->size[0];
    r0->size[0] = bsub - b_b;
    emxEnsureCapacity(sp, (emxArray__common *)r0, nc1, sizeof(int32_T),
                      &emlrtRTEI);
    loop_ub = bsub - b_b;
    for (nc1 = 0; nc1 < loop_ub; nc1++) {
      r0->data[nc1] = b_b + nc1;
    }

    loop_ub = neighbors->size[1];
    b_b = r1->size[0];
    r1->size[0] = loop_ub;
    emxEnsureCapacity(sp, (emxArray__common *)r1, b_b, sizeof(int32_T),
                      &emlrtRTEI);
    for (b_b = 0; b_b < loop_ub; b_b++) {
      r1->data[b_b] = b_b;
    }

    loop_ub = idx->size[1];
    b_b = pts->size[0] * pts->size[1];
    pts->size[0] = loop_ub;
    pts->size[1] = ck - i0;
    emxEnsureCapacity(sp, (emxArray__common *)pts, b_b, sizeof(real_T),
                      &emlrtRTEI);
    csz_idx_1 = ck - i0;
    for (b_b = 0; b_b < csz_idx_1; b_b++) {
      for (ck = 0; ck < loop_ub; ck++) {
        pts->data[ck + pts->size[0] * b_b] = idx->data[(i0 + b_b) + idx->size[0]
          * ck];
      }
    }

    iv0[0] = r0->size[0];
    iv0[1] = r1->size[0];
    emlrtSubAssignSizeCheckR2012b(iv0, 2, *(int32_T (*)[2])pts->size, 2,
      &b_emlrtECI, sp);
    loop_ub = pts->size[1];
    for (i0 = 0; i0 < loop_ub; i0++) {
      csz_idx_1 = pts->size[0];
      for (b_b = 0; b_b < csz_idx_1; b_b++) {
        neighbors->data[r0->data[b_b] + neighbors->size[0] * r1->data[i0]] =
          pts->data[b_b + pts->size[0] * i0];
      }
    }

    if (2.0 > k + 1.0) {
      i0 = 0;
      ck = 0;
    } else {
      i0 = b->size[0];
      if (!(2 <= i0)) {
        emlrtDynamicBoundsCheckR2012b(2, 1, i0, &g_emlrtBCI, sp);
      }

      i0 = 1;
      b_b = b->size[0];
      ck = (int32_T)(k + 1.0);
      if (!((ck >= 1) && (ck <= b_b))) {
        emlrtDynamicBoundsCheckR2012b(ck, 1, b_b, &h_emlrtBCI, sp);
      }
    }

    b_b = start->size[1];
    nc1 = i + 1;
    if (!((nc1 >= 1) && (nc1 <= b_b))) {
      emlrtDynamicBoundsCheckR2012b(nc1, 1, b_b, &m_emlrtBCI, sp);
    }

    b_b = stop->size[1];
    nc1 = i + 1;
    if (!((nc1 >= 1) && (nc1 <= b_b))) {
      emlrtDynamicBoundsCheckR2012b(nc1, 1, b_b, &n_emlrtBCI, sp);
    }

    if (start->data[i] > stop->data[i]) {
      b_b = 0;
      bsub = 0;
    } else {
      alpha1 = start->data[i];
      if (alpha1 != (int32_T)muDoubleScalarFloor(alpha1)) {
        emlrtIntegerCheckR2012b(alpha1, &e_emlrtDCI, sp);
      }

      b_b = distances->size[0];
      nc1 = (int32_T)alpha1;
      if (!((nc1 >= 1) && (nc1 <= b_b))) {
        emlrtDynamicBoundsCheckR2012b(nc1, 1, b_b, &o_emlrtBCI, sp);
      }

      b_b = nc1 - 1;
      alpha1 = stop->data[i];
      if (alpha1 != (int32_T)muDoubleScalarFloor(alpha1)) {
        emlrtIntegerCheckR2012b(alpha1, &f_emlrtDCI, sp);
      }

      nc1 = distances->size[0];
      bsub = (int32_T)alpha1;
      if (!((bsub >= 1) && (bsub <= nc1))) {
        emlrtDynamicBoundsCheckR2012b(bsub, 1, nc1, &p_emlrtBCI, sp);
      }
    }

    nc1 = r0->size[0];
    r0->size[0] = bsub - b_b;
    emxEnsureCapacity(sp, (emxArray__common *)r0, nc1, sizeof(int32_T),
                      &emlrtRTEI);
    loop_ub = bsub - b_b;
    for (nc1 = 0; nc1 < loop_ub; nc1++) {
      r0->data[nc1] = b_b + nc1;
    }

    loop_ub = distances->size[1];
    b_b = r1->size[0];
    r1->size[0] = loop_ub;
    emxEnsureCapacity(sp, (emxArray__common *)r1, b_b, sizeof(int32_T),
                      &emlrtRTEI);
    for (b_b = 0; b_b < loop_ub; b_b++) {
      r1->data[b_b] = b_b;
    }

    loop_ub = b->size[1];
    b_b = pts->size[0] * pts->size[1];
    pts->size[0] = loop_ub;
    pts->size[1] = ck - i0;
    emxEnsureCapacity(sp, (emxArray__common *)pts, b_b, sizeof(real_T),
                      &emlrtRTEI);
    csz_idx_1 = ck - i0;
    for (b_b = 0; b_b < csz_idx_1; b_b++) {
      for (ck = 0; ck < loop_ub; ck++) {
        pts->data[ck + pts->size[0] * b_b] = b->data[(i0 + b_b) + b->size[0] *
          ck];
      }
    }

    iv0[0] = r0->size[0];
    iv0[1] = r1->size[0];
    emlrtSubAssignSizeCheckR2012b(iv0, 2, *(int32_T (*)[2])pts->size, 2,
      &c_emlrtECI, sp);
    loop_ub = pts->size[1];
    for (i0 = 0; i0 < loop_ub; i0++) {
      csz_idx_1 = pts->size[0];
      for (b_b = 0; b_b < csz_idx_1; b_b++) {
        distances->data[r0->data[b_b] + distances->size[0] * r1->data[i0]] =
          pts->data[b_b + pts->size[0] * i0];
      }
    }

    i++;
    if (*emlrtBreakCheckR2012bFlagVar != 0) {
      emlrtBreakCheckR2012b(sp);
    }
  }

  emxFree_real_T(&r2);
  emxFree_int32_T(&iidx);
  emxFree_real_T(&cv);
  emxFree_real_T(&bv);
  emxFree_real_T(&av);
  emxFree_real_T(&a);
  emxFree_real_T(&b);
  emxFree_int32_T(&r1);
  emxFree_int32_T(&r0);
  emxFree_int32_T(&idx);
  emxFree_real_T(&stop);
  emxFree_real_T(&start);
  emxFree_real_T(&S);
  emxFree_real_T(&pts);
  emlrtHeapReferenceStackLeaveFcnR2012b(sp);
}

/* End of code generation (nn_batch.c) */
