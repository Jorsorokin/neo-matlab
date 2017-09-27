/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * find_blockIndex_range.c
 *
 * Code generation for function 'find_blockIndex_range'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "nn_batch.h"
#include "find_blockIndex_range.h"
#include "nn_batch_emxutil.h"

/* Variable Definitions */
static emlrtRSInfo f_emlrtRSI = { 21,  /* lineNo */
  "find_blockIndex_range",             /* fcnName */
  "Z:\\matlabscripts\\SpikeSorting\\Functions\\misc\\find_blockIndex_range.m"/* pathName */
};

static emlrtRSInfo g_emlrtRSI = { 22,  /* lineNo */
  "find_blockIndex_range",             /* fcnName */
  "Z:\\matlabscripts\\SpikeSorting\\Functions\\misc\\find_blockIndex_range.m"/* pathName */
};

static emlrtRSInfo h_emlrtRSI = { 98,  /* lineNo */
  "colon",                             /* fcnName */
  "C:\\MATLAB\\R2017a\\toolbox\\eml\\lib\\matlab\\ops\\colon.m"/* pathName */
};

static emlrtRSInfo i_emlrtRSI = { 282, /* lineNo */
  "colon",                             /* fcnName */
  "C:\\MATLAB\\R2017a\\toolbox\\eml\\lib\\matlab\\ops\\colon.m"/* pathName */
};

static emlrtRSInfo j_emlrtRSI = { 290, /* lineNo */
  "colon",                             /* fcnName */
  "C:\\MATLAB\\R2017a\\toolbox\\eml\\lib\\matlab\\ops\\colon.m"/* pathName */
};

static emlrtRTEInfo g_emlrtRTEI = { 1, /* lineNo */
  25,                                  /* colNo */
  "find_blockIndex_range",             /* fName */
  "Z:\\matlabscripts\\SpikeSorting\\Functions\\misc\\find_blockIndex_range.m"/* pName */
};

static emlrtRTEInfo h_emlrtRTEI = { 21,/* lineNo */
  5,                                   /* colNo */
  "find_blockIndex_range",             /* fName */
  "Z:\\matlabscripts\\SpikeSorting\\Functions\\misc\\find_blockIndex_range.m"/* pName */
};

static emlrtRTEInfo i_emlrtRTEI = { 22,/* lineNo */
  5,                                   /* colNo */
  "find_blockIndex_range",             /* fName */
  "Z:\\matlabscripts\\SpikeSorting\\Functions\\misc\\find_blockIndex_range.m"/* pName */
};

static emlrtRTEInfo x_emlrtRTEI = { 388,/* lineNo */
  15,                                  /* colNo */
  "colon",                             /* fName */
  "C:\\MATLAB\\R2017a\\toolbox\\eml\\lib\\matlab\\ops\\colon.m"/* pName */
};

static emlrtBCInfo q_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  24,                                  /* lineNo */
  38,                                  /* colNo */
  "tempStart",                         /* aName */
  "find_blockIndex_range",             /* fName */
  "Z:\\matlabscripts\\SpikeSorting\\Functions\\misc\\find_blockIndex_range.m",/* pName */
  0                                    /* checkKind */
};

/* Function Definitions */
void find_blockIndex_range(const emlrtStack *sp, real_T n, real_T m,
  emxArray_real_T *start, emxArray_real_T *stop)
{
  real_T maxPts;
  int32_T k;
  real_T b_remainder;
  real_T kd;
  emxArray_real_T *tempStart;
  real_T ndbl;
  real_T apnd;
  real_T cdiff;
  int32_T nm1d2;
  emxArray_real_T *tempStop;
  real_T absb;
  int32_T b_n;
  real_T absa;
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

  /*  [start,stop] = find_indexing_range( n,m,(maxSize) ) */
  /*  */
  /*  find the indexing ranges to perform block looping over points "1:m" for  */
  /*  nearest-neighbor, etc. This greatly improves performance (not looping over  */
  /*  all points) while avoiding massive (> 1 GB) matrices. Default "maxSize", */
  /*  which refers to the maximum matrix size allowed, is 100,000,000 entries. */
  /*  */
  /*  "n" refers to number of elements in first vector, "m" refers to number of */
  /*  elements in second. This is more flexible than simply providing one */
  /*  number, as neighborhood graphs etc. may not necessarily be square. */
  /*  get max Pts possible */
  maxPts = muDoubleScalarFloor(1.0E+8 / n);
  if (maxPts < m) {
    b_remainder = m;
    if (!muDoubleScalarIsInf(maxPts)) {
      if (m == 0.0) {
        b_remainder = maxPts * 0.0;
      } else {
        if (maxPts != 0.0) {
          b_remainder = muDoubleScalarRem(m, maxPts);
          if (b_remainder == 0.0) {
            b_remainder = maxPts * 0.0;
          } else {
            if (maxPts < 0.0) {
              b_remainder += maxPts;
            }
          }
        }
      }
    } else {
      if (maxPts != 0.0) {
        b_remainder = rtNaN;
      }
    }

    st.site = &f_emlrtRSI;
    kd = m - b_remainder;
    emxInit_real_T(&st, &tempStart, 2, &h_emlrtRTEI, true);
    if (muDoubleScalarIsNaN(kd)) {
      k = tempStart->size[0] * tempStart->size[1];
      tempStart->size[0] = 1;
      tempStart->size[1] = 1;
      emxEnsureCapacity(&st, (emxArray__common *)tempStart, k, sizeof(real_T),
                        &g_emlrtRTEI);
      tempStart->data[0] = rtNaN;
    } else if ((maxPts == 0.0) || ((1.0 < kd) && (maxPts < 0.0)) || ((kd < 1.0) &&
                (maxPts > 0.0))) {
      k = tempStart->size[0] * tempStart->size[1];
      tempStart->size[0] = 1;
      tempStart->size[1] = 0;
      emxEnsureCapacity(&st, (emxArray__common *)tempStart, k, sizeof(real_T),
                        &g_emlrtRTEI);
    } else if (muDoubleScalarIsInf(kd) && (muDoubleScalarIsInf(maxPts) || (1.0 ==
      kd))) {
      k = tempStart->size[0] * tempStart->size[1];
      tempStart->size[0] = 1;
      tempStart->size[1] = 1;
      emxEnsureCapacity(&st, (emxArray__common *)tempStart, k, sizeof(real_T),
                        &g_emlrtRTEI);
      tempStart->data[0] = rtNaN;
    } else if (muDoubleScalarIsInf(maxPts)) {
      k = tempStart->size[0] * tempStart->size[1];
      tempStart->size[0] = 1;
      tempStart->size[1] = 1;
      emxEnsureCapacity(&st, (emxArray__common *)tempStart, k, sizeof(real_T),
                        &g_emlrtRTEI);
      tempStart->data[0] = 1.0;
    } else if (maxPts == maxPts) {
      k = tempStart->size[0] * tempStart->size[1];
      tempStart->size[0] = 1;
      tempStart->size[1] = (int32_T)muDoubleScalarFloor((kd - 1.0) / maxPts) + 1;
      emxEnsureCapacity(&st, (emxArray__common *)tempStart, k, sizeof(real_T),
                        &g_emlrtRTEI);
      nm1d2 = (int32_T)muDoubleScalarFloor((kd - 1.0) / maxPts);
      for (k = 0; k <= nm1d2; k++) {
        tempStart->data[tempStart->size[0] * k] = 1.0 + maxPts * (real_T)k;
      }
    } else {
      b_st.site = &h_emlrtRSI;
      ndbl = muDoubleScalarFloor((kd - 1.0) / maxPts + 0.5);
      apnd = 1.0 + ndbl * maxPts;
      if (maxPts > 0.0) {
        cdiff = apnd - kd;
      } else {
        cdiff = kd - apnd;
      }

      absb = muDoubleScalarAbs(kd);
      if (muDoubleScalarAbs(cdiff) < 4.4408920985006262E-16 * muDoubleScalarMax
          (1.0, absb)) {
        ndbl++;
        apnd = kd;
      } else if (cdiff > 0.0) {
        apnd = 1.0 + (ndbl - 1.0) * maxPts;
      } else {
        ndbl++;
      }

      if (ndbl >= 0.0) {
        b_n = (int32_T)ndbl;
      } else {
        b_n = 0;
      }

      c_st.site = &i_emlrtRSI;
      if (2.147483647E+9 < ndbl) {
        emlrtErrorWithMessageIdR2012b(&c_st, &x_emlrtRTEI,
          "Coder:MATLAB:pmaxsize", 0);
      }

      k = tempStart->size[0] * tempStart->size[1];
      tempStart->size[0] = 1;
      tempStart->size[1] = b_n;
      emxEnsureCapacity(&b_st, (emxArray__common *)tempStart, k, sizeof(real_T),
                        &g_emlrtRTEI);
      if (b_n > 0) {
        tempStart->data[0] = 1.0;
        if (b_n > 1) {
          tempStart->data[b_n - 1] = apnd;
          nm1d2 = (b_n - 1) / 2;
          c_st.site = &j_emlrtRSI;
          for (k = 1; k < nm1d2; k++) {
            kd = (real_T)k * maxPts;
            tempStart->data[k] = 1.0 + kd;
            tempStart->data[(b_n - k) - 1] = apnd - kd;
          }

          if (nm1d2 << 1 == b_n - 1) {
            tempStart->data[nm1d2] = (1.0 + apnd) / 2.0;
          } else {
            kd = (real_T)nm1d2 * maxPts;
            tempStart->data[nm1d2] = 1.0 + kd;
            tempStart->data[nm1d2 + 1] = apnd - kd;
          }
        }
      }
    }

    st.site = &g_emlrtRSI;
    kd = m - b_remainder;
    emxInit_real_T(&st, &tempStop, 2, &i_emlrtRTEI, true);
    if (muDoubleScalarIsNaN(kd)) {
      k = tempStop->size[0] * tempStop->size[1];
      tempStop->size[0] = 1;
      tempStop->size[1] = 1;
      emxEnsureCapacity(&st, (emxArray__common *)tempStop, k, sizeof(real_T),
                        &g_emlrtRTEI);
      tempStop->data[0] = rtNaN;
    } else if ((maxPts == 0.0) || ((maxPts < kd) && (maxPts < 0.0)) || ((kd <
                 maxPts) && (maxPts > 0.0))) {
      k = tempStop->size[0] * tempStop->size[1];
      tempStop->size[0] = 1;
      tempStop->size[1] = 0;
      emxEnsureCapacity(&st, (emxArray__common *)tempStop, k, sizeof(real_T),
                        &g_emlrtRTEI);
    } else if ((muDoubleScalarIsInf(maxPts) || muDoubleScalarIsInf(kd)) &&
               (muDoubleScalarIsInf(maxPts) || (maxPts == kd))) {
      k = tempStop->size[0] * tempStop->size[1];
      tempStop->size[0] = 1;
      tempStop->size[1] = 1;
      emxEnsureCapacity(&st, (emxArray__common *)tempStop, k, sizeof(real_T),
                        &g_emlrtRTEI);
      tempStop->data[0] = rtNaN;
    } else if (muDoubleScalarIsInf(maxPts)) {
      k = tempStop->size[0] * tempStop->size[1];
      tempStop->size[0] = 1;
      tempStop->size[1] = 1;
      emxEnsureCapacity(&st, (emxArray__common *)tempStop, k, sizeof(real_T),
                        &g_emlrtRTEI);
      tempStop->data[0] = maxPts;
    } else if (maxPts == maxPts) {
      k = tempStop->size[0] * tempStop->size[1];
      tempStop->size[0] = 1;
      tempStop->size[1] = (int32_T)muDoubleScalarFloor((kd - maxPts) / maxPts) +
        1;
      emxEnsureCapacity(&st, (emxArray__common *)tempStop, k, sizeof(real_T),
                        &g_emlrtRTEI);
      nm1d2 = (int32_T)muDoubleScalarFloor((kd - maxPts) / maxPts);
      for (k = 0; k <= nm1d2; k++) {
        tempStop->data[tempStop->size[0] * k] = maxPts + maxPts * (real_T)k;
      }
    } else {
      b_st.site = &h_emlrtRSI;
      ndbl = muDoubleScalarFloor((kd - maxPts) / maxPts + 0.5);
      apnd = maxPts + ndbl * maxPts;
      if (maxPts > 0.0) {
        cdiff = apnd - kd;
      } else {
        cdiff = kd - apnd;
      }

      absa = muDoubleScalarAbs(maxPts);
      absb = muDoubleScalarAbs(kd);
      if (muDoubleScalarAbs(cdiff) < 4.4408920985006262E-16 * muDoubleScalarMax
          (absa, absb)) {
        ndbl++;
        apnd = kd;
      } else if (cdiff > 0.0) {
        apnd = maxPts + (ndbl - 1.0) * maxPts;
      } else {
        ndbl++;
      }

      if (ndbl >= 0.0) {
        b_n = (int32_T)ndbl;
      } else {
        b_n = 0;
      }

      c_st.site = &i_emlrtRSI;
      if (2.147483647E+9 < ndbl) {
        emlrtErrorWithMessageIdR2012b(&c_st, &x_emlrtRTEI,
          "Coder:MATLAB:pmaxsize", 0);
      }

      k = tempStop->size[0] * tempStop->size[1];
      tempStop->size[0] = 1;
      tempStop->size[1] = b_n;
      emxEnsureCapacity(&b_st, (emxArray__common *)tempStop, k, sizeof(real_T),
                        &g_emlrtRTEI);
      if (b_n > 0) {
        tempStop->data[0] = maxPts;
        if (b_n > 1) {
          tempStop->data[b_n - 1] = apnd;
          nm1d2 = (b_n - 1) / 2;
          c_st.site = &j_emlrtRSI;
          for (k = 1; k < nm1d2; k++) {
            kd = (real_T)k * maxPts;
            tempStop->data[k] = maxPts + kd;
            tempStop->data[(b_n - k) - 1] = apnd - kd;
          }

          if (nm1d2 << 1 == b_n - 1) {
            tempStop->data[nm1d2] = (maxPts + apnd) / 2.0;
          } else {
            kd = (real_T)nm1d2 * maxPts;
            tempStop->data[nm1d2] = maxPts + kd;
            tempStop->data[nm1d2 + 1] = apnd - kd;
          }
        }
      }
    }

    if (b_remainder > 0.0) {
      k = start->size[0] * start->size[1];
      start->size[0] = 1;
      start->size[1] = tempStart->size[1] + 1;
      emxEnsureCapacity(sp, (emxArray__common *)start, k, sizeof(real_T),
                        &g_emlrtRTEI);
      nm1d2 = tempStart->size[1];
      for (k = 0; k < nm1d2; k++) {
        start->data[start->size[0] * k] = tempStart->data[tempStart->size[0] * k];
      }

      k = tempStart->size[1];
      nm1d2 = tempStart->size[1];
      if (!((nm1d2 >= 1) && (nm1d2 <= k))) {
        emlrtDynamicBoundsCheckR2012b(nm1d2, 1, k, &q_emlrtBCI, sp);
      }

      start->data[start->size[0] * tempStart->size[1]] = tempStart->data[nm1d2 -
        1] + b_remainder;
      k = stop->size[0] * stop->size[1];
      stop->size[0] = 1;
      stop->size[1] = tempStop->size[1] + 1;
      emxEnsureCapacity(sp, (emxArray__common *)stop, k, sizeof(real_T),
                        &g_emlrtRTEI);
      nm1d2 = tempStop->size[1];
      for (k = 0; k < nm1d2; k++) {
        stop->data[stop->size[0] * k] = tempStop->data[tempStop->size[0] * k];
      }

      stop->data[stop->size[0] * tempStop->size[1]] = m;
    } else {
      k = start->size[0] * start->size[1];
      start->size[0] = 1;
      start->size[1] = tempStart->size[1];
      emxEnsureCapacity(sp, (emxArray__common *)start, k, sizeof(real_T),
                        &g_emlrtRTEI);
      nm1d2 = tempStart->size[0] * tempStart->size[1];
      for (k = 0; k < nm1d2; k++) {
        start->data[k] = tempStart->data[k];
      }

      k = stop->size[0] * stop->size[1];
      stop->size[0] = 1;
      stop->size[1] = tempStop->size[1];
      emxEnsureCapacity(sp, (emxArray__common *)stop, k, sizeof(real_T),
                        &g_emlrtRTEI);
      nm1d2 = tempStop->size[0] * tempStop->size[1];
      for (k = 0; k < nm1d2; k++) {
        stop->data[k] = tempStop->data[k];
      }
    }

    emxFree_real_T(&tempStop);
    emxFree_real_T(&tempStart);
  } else {
    k = start->size[0] * start->size[1];
    start->size[0] = 1;
    start->size[1] = 1;
    emxEnsureCapacity(sp, (emxArray__common *)start, k, sizeof(real_T),
                      &g_emlrtRTEI);
    start->data[0] = 1.0;
    k = stop->size[0] * stop->size[1];
    stop->size[0] = 1;
    stop->size[1] = 1;
    emxEnsureCapacity(sp, (emxArray__common *)stop, k, sizeof(real_T),
                      &g_emlrtRTEI);
    stop->data[0] = m;
  }

  emlrtHeapReferenceStackLeaveFcnR2012b(sp);
}

/* End of code generation (find_blockIndex_range.c) */
