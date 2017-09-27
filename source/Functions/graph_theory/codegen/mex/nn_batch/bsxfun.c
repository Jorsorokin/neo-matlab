/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * bsxfun.c
 *
 * Code generation for function 'bsxfun'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "nn_batch.h"
#include "bsxfun.h"
#include "nn_batch_emxutil.h"
#include "eml_int_forloop_overflow_check.h"
#include "nn_batch_data.h"

/* Variable Definitions */
static emlrtRTEInfo k_emlrtRTEI = { 1, /* lineNo */
  14,                                  /* colNo */
  "bsxfun",                            /* fName */
  "C:\\MATLAB\\R2017a\\toolbox\\eml\\lib\\matlab\\elmat\\bsxfun.m"/* pName */
};

/* Function Definitions */
void bsxfun(const emlrtStack *sp, const emxArray_real_T *a, const
            emxArray_real_T *b, emxArray_real_T *c)
{
  int32_T na1;
  int32_T nb1;
  boolean_T overflow;
  int32_T b_c;
  int32_T csz_idx_1;
  int32_T i1;
  emxArray_real_T *av;
  emxArray_real_T *bv;
  uint32_T a_idx_0;
  int32_T bk;
  int32_T nc1;
  int32_T b_b;
  int32_T ck;
  emxArray_real_T *cv;
  emlrtStack st;
  emlrtStack b_st;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  emlrtHeapReferenceStackEnterFcnR2012b(sp);
  na1 = a->size[0];
  nb1 = b->size[0];
  if ((a->size[0] != 1) && (b->size[0] != 1) && (a->size[0] != b->size[0])) {
    overflow = false;
  } else {
    overflow = true;
  }

  if (!overflow) {
    emlrtErrorWithMessageIdR2012b(sp, &w_emlrtRTEI,
      "MATLAB:bsxfun:arrayDimensionsMustMatch", 0);
  }

  overflow = !(a->size[0] != b->size[0]);
  if (!overflow) {
    emlrtErrorWithMessageIdR2012b(sp, &v_emlrtRTEI,
      "Coder:toolbox:bsxfun_dynamicExpansion", 0);
  }

  if (a->size[0] <= b->size[0]) {
    b_c = a->size[0];
  } else {
    b_c = b->size[0];
  }

  csz_idx_1 = b->size[1];
  i1 = c->size[0] * c->size[1];
  c->size[0] = b_c;
  c->size[1] = csz_idx_1;
  emxEnsureCapacity(sp, (emxArray__common *)c, i1, sizeof(real_T), &k_emlrtRTEI);
  if (!((c->size[0] == 0) || (c->size[1] == 0))) {
    emxInit_real_T1(sp, &av, 1, &d_emlrtRTEI, true);
    emxInit_real_T1(sp, &bv, 1, &e_emlrtRTEI, true);
    a_idx_0 = (uint32_T)a->size[0];
    i1 = av->size[0];
    av->size[0] = (int32_T)a_idx_0;
    emxEnsureCapacity(sp, (emxArray__common *)av, i1, sizeof(real_T),
                      &k_emlrtRTEI);
    a_idx_0 = (uint32_T)b->size[0];
    i1 = bv->size[0];
    bv->size[0] = (int32_T)a_idx_0;
    emxEnsureCapacity(sp, (emxArray__common *)bv, i1, sizeof(real_T),
                      &k_emlrtRTEI);
    csz_idx_1 = 1;
    bk = 0;
    nc1 = c->size[0];
    b_b = c->size[0] * c->size[1] - c->size[0];
    st.site = &fb_emlrtRSI;
    overflow = ((!(0 > b_b)) && (b_b > MAX_int32_T - c->size[0]));
    if (overflow) {
      b_st.site = &k_emlrtRSI;
      check_forloop_overflow_error(&b_st);
    }

    ck = 0;
    emxInit_real_T1(sp, &cv, 1, &f_emlrtRTEI, true);
    while (ck <= b_b) {
      st.site = &eb_emlrtRSI;
      if ((!(1 > na1)) && (na1 > 2147483646)) {
        b_st.site = &k_emlrtRSI;
        check_forloop_overflow_error(&b_st);
      }

      for (b_c = 0; b_c + 1 <= na1; b_c++) {
        av->data[b_c] = a->data[b_c];
      }

      st.site = &db_emlrtRSI;
      if ((!(1 > nb1)) && (nb1 > 2147483646)) {
        b_st.site = &k_emlrtRSI;
        check_forloop_overflow_error(&b_st);
      }

      for (b_c = 0; b_c + 1 <= nb1; b_c++) {
        bv->data[b_c] = b->data[bk + b_c];
      }

      i1 = cv->size[0];
      cv->size[0] = av->size[0];
      emxEnsureCapacity(sp, (emxArray__common *)cv, i1, sizeof(real_T),
                        &k_emlrtRTEI);
      b_c = av->size[0];
      for (i1 = 0; i1 < b_c; i1++) {
        cv->data[i1] = av->data[i1] - bv->data[i1];
      }

      st.site = &cb_emlrtRSI;
      if (nc1 > 2147483646) {
        b_st.site = &k_emlrtRSI;
        check_forloop_overflow_error(&b_st);
      }

      for (b_c = 0; b_c + 1 <= nc1; b_c++) {
        c->data[ck + b_c] = cv->data[b_c];
      }

      if (csz_idx_1 < b->size[1]) {
        bk += nb1;
        csz_idx_1++;
      } else {
        csz_idx_1 = 1;
      }

      ck += nc1;
    }

    emxFree_real_T(&cv);
    emxFree_real_T(&bv);
    emxFree_real_T(&av);
  }

  emlrtHeapReferenceStackLeaveFcnR2012b(sp);
}

/* End of code generation (bsxfun.c) */
