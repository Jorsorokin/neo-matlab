/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * dot.c
 *
 * Code generation for function 'dot'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "nn_batch.h"
#include "dot.h"
#include "nn_batch_emxutil.h"
#include "eml_int_forloop_overflow_check.h"
#include "nn_batch_data.h"
#include "blas.h"

/* Variable Definitions */
static emlrtRSInfo l_emlrtRSI = { 50,  /* lineNo */
  "dot",                               /* fcnName */
  "C:\\MATLAB\\R2017a\\toolbox\\eml\\lib\\matlab\\specfun\\dot.m"/* pathName */
};

static emlrtRSInfo m_emlrtRSI = { 67,  /* lineNo */
  "dot",                               /* fcnName */
  "C:\\MATLAB\\R2017a\\toolbox\\eml\\lib\\matlab\\specfun\\dot.m"/* pathName */
};

static emlrtRSInfo n_emlrtRSI = { 71,  /* lineNo */
  "dot",                               /* fcnName */
  "C:\\MATLAB\\R2017a\\toolbox\\eml\\lib\\matlab\\specfun\\dot.m"/* pathName */
};

static emlrtRTEInfo j_emlrtRTEI = { 1, /* lineNo */
  14,                                  /* colNo */
  "dot",                               /* fName */
  "C:\\MATLAB\\R2017a\\toolbox\\eml\\lib\\matlab\\specfun\\dot.m"/* pName */
};

static emlrtRTEInfo ab_emlrtRTEI = { 13,/* lineNo */
  15,                                  /* colNo */
  "dot",                               /* fName */
  "C:\\MATLAB\\R2017a\\toolbox\\eml\\lib\\matlab\\specfun\\dot.m"/* pName */
};

/* Function Definitions */
void dot(const emlrtStack *sp, const emxArray_real_T *a, const emxArray_real_T
         *b, emxArray_real_T *c)
{
  int32_T vlen;
  uint32_T varargin_1[2];
  boolean_T overflow;
  uint32_T varargin_2[2];
  boolean_T p;
  boolean_T exitg1;
  int32_T vstride;
  int32_T ic;
  int32_T i1;
  int32_T j;
  real_T b_c;
  ptrdiff_t n_t;
  ptrdiff_t incx_t;
  ptrdiff_t incy_t;
  emlrtStack st;
  emlrtStack b_st;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  for (vlen = 0; vlen < 2; vlen++) {
    varargin_1[vlen] = (uint32_T)a->size[vlen];
  }

  for (vlen = 0; vlen < 2; vlen++) {
    varargin_2[vlen] = (uint32_T)b->size[vlen];
  }

  overflow = false;
  p = true;
  vlen = 0;
  exitg1 = false;
  while ((!exitg1) && (vlen < 2)) {
    if (!((int32_T)varargin_1[vlen] == (int32_T)varargin_2[vlen])) {
      p = false;
      exitg1 = true;
    } else {
      vlen++;
    }
  }

  if (p) {
    overflow = true;
  }

  if (overflow || (((a->size[0] == 1) || (a->size[1] == 1)) && ((b->size[0] == 1)
        || (b->size[1] == 1)) && (a->size[0] * a->size[1] == b->size[0] *
        b->size[1]))) {
    overflow = true;
  } else {
    overflow = false;
  }

  if (!overflow) {
    emlrtErrorWithMessageIdR2012b(sp, &ab_emlrtRTEI,
      "MATLAB:dot:InputSizeMismatch", 0);
  }

  if (a->size[1] == 1) {
    vlen = c->size[0];
    c->size[0] = a->size[0];
    emxEnsureCapacity(sp, (emxArray__common *)c, vlen, sizeof(real_T),
                      &j_emlrtRTEI);
    vlen = c->size[0];
    st.site = &l_emlrtRSI;
    overflow = ((!(1 > c->size[0])) && (c->size[0] > 2147483646));
    if (overflow) {
      b_st.site = &k_emlrtRSI;
      check_forloop_overflow_error(&b_st);
    }

    for (j = 0; j + 1 <= vlen; j++) {
      c->data[j] = a->data[j] * b->data[j];
    }
  } else {
    vlen = c->size[0];
    c->size[0] = a->size[0];
    emxEnsureCapacity(sp, (emxArray__common *)c, vlen, sizeof(real_T),
                      &j_emlrtRTEI);
    vlen = a->size[1];
    vstride = a->size[0];
    ic = -1;
    i1 = -1;
    st.site = &m_emlrtRSI;
    overflow = ((!(1 > a->size[0])) && (a->size[0] > 2147483646));
    if (overflow) {
      b_st.site = &k_emlrtRSI;
      check_forloop_overflow_error(&b_st);
    }

    for (j = 1; j <= vstride; j++) {
      ic++;
      i1++;
      st.site = &n_emlrtRSI;
      if (vlen < 1) {
        b_c = 0.0;
      } else {
        n_t = (ptrdiff_t)vlen;
        incx_t = (ptrdiff_t)vstride;
        incy_t = (ptrdiff_t)vstride;
        b_c = ddot(&n_t, &a->data[i1], &incx_t, &b->data[i1], &incy_t);
      }

      c->data[ic] = b_c;
    }
  }
}

/* End of code generation (dot.c) */
