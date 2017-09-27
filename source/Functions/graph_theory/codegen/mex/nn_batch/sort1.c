/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * sort1.c
 *
 * Code generation for function 'sort1'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "nn_batch.h"
#include "sort1.h"
#include "nn_batch_emxutil.h"
#include "eml_int_forloop_overflow_check.h"
#include "sortIdx.h"
#include "nn_batch_data.h"

/* Variable Definitions */
static emlrtRSInfo lb_emlrtRSI = { 19, /* lineNo */
  "sort",                              /* fcnName */
  "C:\\MATLAB\\R2017a\\toolbox\\eml\\eml\\+coder\\+internal\\sort.m"/* pathName */
};

static emlrtRSInfo mb_emlrtRSI = { 84, /* lineNo */
  "sort",                              /* fcnName */
  "C:\\MATLAB\\R2017a\\toolbox\\eml\\eml\\+coder\\+internal\\sort.m"/* pathName */
};

static emlrtRSInfo nb_emlrtRSI = { 82, /* lineNo */
  "sort",                              /* fcnName */
  "C:\\MATLAB\\R2017a\\toolbox\\eml\\eml\\+coder\\+internal\\sort.m"/* pathName */
};

static emlrtRSInfo ob_emlrtRSI = { 79, /* lineNo */
  "sort",                              /* fcnName */
  "C:\\MATLAB\\R2017a\\toolbox\\eml\\eml\\+coder\\+internal\\sort.m"/* pathName */
};

static emlrtRSInfo pb_emlrtRSI = { 76, /* lineNo */
  "sort",                              /* fcnName */
  "C:\\MATLAB\\R2017a\\toolbox\\eml\\eml\\+coder\\+internal\\sort.m"/* pathName */
};

static emlrtRSInfo qb_emlrtRSI = { 74, /* lineNo */
  "sort",                              /* fcnName */
  "C:\\MATLAB\\R2017a\\toolbox\\eml\\eml\\+coder\\+internal\\sort.m"/* pathName */
};

static emlrtRSInfo rb_emlrtRSI = { 71, /* lineNo */
  "sort",                              /* fcnName */
  "C:\\MATLAB\\R2017a\\toolbox\\eml\\eml\\+coder\\+internal\\sort.m"/* pathName */
};

static emlrtRSInfo sb_emlrtRSI = { 50, /* lineNo */
  "sort",                              /* fcnName */
  "C:\\MATLAB\\R2017a\\toolbox\\eml\\eml\\+coder\\+internal\\sort.m"/* pathName */
};

static emlrtRSInfo tb_emlrtRSI = { 49, /* lineNo */
  "prodsize",                          /* fcnName */
  "C:\\MATLAB\\R2017a\\toolbox\\eml\\eml\\+coder\\+internal\\prodsize.m"/* pathName */
};

static emlrtRTEInfo l_emlrtRTEI = { 1, /* lineNo */
  20,                                  /* colNo */
  "sort",                              /* fName */
  "C:\\MATLAB\\R2017a\\toolbox\\eml\\eml\\+coder\\+internal\\sort.m"/* pName */
};

static emlrtRTEInfo o_emlrtRTEI = { 52,/* lineNo */
  1,                                   /* colNo */
  "sort",                              /* fName */
  "C:\\MATLAB\\R2017a\\toolbox\\eml\\eml\\+coder\\+internal\\sort.m"/* pName */
};

/* Function Declarations */
static void b_sort(const emlrtStack *sp, emxArray_real_T *x, int32_T dim,
                   emxArray_int32_T *idx);

/* Function Definitions */
static void b_sort(const emlrtStack *sp, emxArray_real_T *x, int32_T dim,
                   emxArray_int32_T *idx)
{
  emxArray_real_T *vwork;
  int32_T i3;
  int32_T vstride;
  int32_T k;
  int32_T iv4[2];
  int32_T npages;
  int32_T pagesize;
  int32_T i;
  emxArray_int32_T *iidx;
  int32_T pageoffset;
  int32_T j;
  int32_T idx0;
  emlrtStack st;
  emlrtStack b_st;
  emlrtStack c_st;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  emlrtHeapReferenceStackEnterFcnR2012b(sp);
  emxInit_real_T1(sp, &vwork, 1, &o_emlrtRTEI, true);
  st.site = &sb_emlrtRSI;
  i3 = x->size[dim - 1];
  vstride = x->size[dim - 1];
  k = vwork->size[0];
  vwork->size[0] = vstride;
  emxEnsureCapacity(sp, (emxArray__common *)vwork, k, sizeof(real_T),
                    &l_emlrtRTEI);
  for (k = 0; k < 2; k++) {
    iv4[k] = x->size[k];
  }

  k = idx->size[0] * idx->size[1];
  idx->size[0] = iv4[0];
  idx->size[1] = iv4[1];
  emxEnsureCapacity(sp, (emxArray__common *)idx, k, sizeof(int32_T),
                    &l_emlrtRTEI);
  st.site = &rb_emlrtRSI;
  vstride = 1;
  b_st.site = &tb_emlrtRSI;
  if ((!(1 > dim - 1)) && (dim - 1 > 2147483646)) {
    c_st.site = &k_emlrtRSI;
    check_forloop_overflow_error(&c_st);
  }

  k = 1;
  while (k <= dim - 1) {
    vstride *= x->size[0];
    k = 2;
  }

  npages = 1;
  k = dim + 1;
  while (k < 3) {
    npages *= x->size[1];
    k = 3;
  }

  pagesize = x->size[dim - 1] * vstride;
  st.site = &qb_emlrtRSI;
  if ((!(1 > npages)) && (npages > 2147483646)) {
    b_st.site = &k_emlrtRSI;
    check_forloop_overflow_error(&b_st);
  }

  i = 1;
  emxInit_int32_T1(sp, &iidx, 1, &l_emlrtRTEI, true);
  while (i <= npages) {
    pageoffset = (i - 1) * pagesize;
    st.site = &pb_emlrtRSI;
    if ((!(1 > vstride)) && (vstride > 2147483646)) {
      b_st.site = &k_emlrtRSI;
      check_forloop_overflow_error(&b_st);
    }

    for (j = 0; j + 1 <= vstride; j++) {
      idx0 = pageoffset + j;
      st.site = &ob_emlrtRSI;
      if ((!(1 > i3)) && (i3 > 2147483646)) {
        b_st.site = &k_emlrtRSI;
        check_forloop_overflow_error(&b_st);
      }

      for (k = 0; k + 1 <= i3; k++) {
        vwork->data[k] = x->data[idx0 + k * vstride];
      }

      st.site = &nb_emlrtRSI;
      sortIdx(&st, vwork, iidx);
      st.site = &mb_emlrtRSI;
      for (k = 0; k + 1 <= i3; k++) {
        x->data[idx0 + k * vstride] = vwork->data[k];
        idx->data[idx0 + k * vstride] = iidx->data[k];
      }
    }

    i++;
  }

  emxFree_int32_T(&iidx);
  emxFree_real_T(&vwork);
  emlrtHeapReferenceStackLeaveFcnR2012b(sp);
}

void sort(const emlrtStack *sp, emxArray_real_T *x, emxArray_int32_T *idx)
{
  int32_T i2;
  emlrtStack st;
  st.prev = sp;
  st.tls = sp->tls;
  i2 = 2;
  if (x->size[0] != 1) {
    i2 = 1;
  }

  st.site = &lb_emlrtRSI;
  b_sort(&st, x, i2, idx);
}

/* End of code generation (sort1.c) */
