/*     CalculiX - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998-2024 Guido Dhondt                     */
/*                                                                        */
/*     This program is free software; you can redistribute it and/or     */
/*     modify it under the terms of the GNU General Public License as    */
/*     published by the Free Software Foundation(version 2);             */
/*                                                                        */
/*     This program is distributed in the hope that it will be useful,   */
/*     but WITHOUT ANY WARRANTY; without even the implied warranty of    */
/*     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the      */
/*     GNU General Public License for more details.                      */
/*                                                                        */
/*     You should have received a copy of the GNU General Public License */
/*     along with this program; if not, write to the Free Software       */
/*     Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.         */

/*
 * Fast sort routines replacing isortii.f, isortid.f, dsort.f.
 *
 * isortii / isortid: LSD radix sort on int64 keys, 8-bit digits (8 passes).
 *   - Signed int64 values are converted to unsigned sort order via sign-bit
 *     XOR before sorting and converted back afterwards.
 *   - Passes where all digit values are identical (cnt[0]==n, meaning the
 *     digit is uniformly 0) are skipped.  FEM node/element indices are
 *     typically < 2^24, so this eliminates 5 of the 8 passes in practice.
 *   - Time: O(8n) worst-case; space: O(n) for two temporary buffers.
 *
 * dsort: struct-packing sort on double keys using the system qsort.
 *   - Key-only path: qsort directly on the key array.
 *   - Key+payload path: pack {double, ITG} pairs, qsort, unpack.
 *   - Time: O(n log n); space: O(n) for the struct array.
 *
 * No external dependencies.  Plain C99.
 *
 * Fortran calling convention (gfortran):
 *   - symbol = lowercase subroutine name + trailing underscore
 *   - all arguments passed by pointer
 *   - integers are int64_t throughout (-fdefault-integer-8)
 *
 * kflag semantics (identical to original SLATEC routines):
 *    1  sort key ascending,  ignore payload
 *   -1  sort key descending, ignore payload
 *    2  sort key ascending,  carry payload
 *   -2  sort key descending, carry payload
 */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/* Fortran integer: always 8 bytes with -fdefault-integer-8 */
typedef long long ITG;

/* ------------------------------------------------------------------ */
/*  Helpers                                                            */
/* ------------------------------------------------------------------ */

/* Map signed int64 to unsigned sort order (flip sign bit).
 * This preserves total order: a < b  <=>  to_u64(a) < to_u64(b). */
static unsigned long long to_u64(ITG x)
{
    return (unsigned long long)x ^ 0x8000000000000000ULL;
}
static ITG from_u64(unsigned long long u)
{
    return (ITG)(u ^ 0x8000000000000000ULL);
}

static void reverse_i(ITG *a, ITG n)
{
    ITG lo = 0, hi = n - 1, tmp;
    while (lo < hi) { tmp = a[lo]; a[lo] = a[hi]; a[hi] = tmp; lo++; hi--; }
}
static void reverse_d(double *a, ITG n)
{
    ITG lo = 0, hi = n - 1;
    double tmp;
    while (lo < hi) { tmp = a[lo]; a[lo] = a[hi]; a[hi] = tmp; lo++; hi--; }
}

/* ------------------------------------------------------------------ */
/*  LSD radix sort: int64 key, optional int64 or double payload        */
/*                                                                     */
/*  Exactly one of vi/vd may be non-NULL (the payload to carry along). */
/* ------------------------------------------------------------------ */
static void lsd_sort_i64(ITG *key, ITG *vi, double *vd, ITG n)
{
    ITG i;
    unsigned long long *uk[2];
    ITG    *ki[2];
    double *kd[2];

    uk[0] = (unsigned long long *)malloc(n * sizeof(unsigned long long));
    uk[1] = (unsigned long long *)malloc(n * sizeof(unsigned long long));
    ki[0] = ki[1] = NULL;
    kd[0] = kd[1] = NULL;

    if (vi) {
        ki[0] = (ITG *)   malloc(n * sizeof(ITG));
        ki[1] = (ITG *)   malloc(n * sizeof(ITG));
    }
    if (vd) {
        kd[0] = (double *)malloc(n * sizeof(double));
        kd[1] = (double *)malloc(n * sizeof(double));
    }

    /* Load into buffer 0 */
    for (i = 0; i < n; i++) uk[0][i] = to_u64(key[i]);
    if (vi) memcpy(ki[0], vi, n * sizeof(ITG));
    if (vd) memcpy(kd[0], vd, n * sizeof(double));

    int cur = 0; /* result lives in buffer cur */

    int pass;
    for (pass = 0; pass < 8; pass++) {
        int shift = pass * 8;
        int nxt   = 1 - cur;
        ITG cnt[256];
        memset(cnt, 0, sizeof(cnt));

        /* Count frequencies of the current 8-bit digit */
        for (i = 0; i < n; i++)
            cnt[(uk[cur][i] >> shift) & 0xFF]++;

        /* If every element has digit value 0, the scatter is a no-op:
         * skip both the scatter and the buffer swap. */
        if (cnt[0] == n) continue;

        /* Exclusive prefix sum */
        ITG total = 0, c;
        for (c = 0; c < 256; c++) {
            ITG old = cnt[c];
            cnt[c]  = total;
            total  += old;
        }

        /* Stable scatter src → dst */
        for (i = 0; i < n; i++) {
            int    bucket = (uk[cur][i] >> shift) & 0xFF;
            ITG    dst    = cnt[bucket]++;
            uk[nxt][dst]  = uk[cur][i];
            if (ki[0]) ki[nxt][dst] = ki[cur][i];
            if (kd[0]) kd[nxt][dst] = kd[cur][i];
        }

        cur = nxt;
    }

    /* Convert back and write results to caller's arrays */
    for (i = 0; i < n; i++) key[i] = from_u64(uk[cur][i]);
    if (vi) memcpy(vi, ki[cur], n * sizeof(ITG));
    if (vd) memcpy(vd, kd[cur], n * sizeof(double));

    free(uk[0]); free(uk[1]);
    if (ki[0]) { free(ki[0]); free(ki[1]); }
    if (kd[0]) { free(kd[0]); free(kd[1]); }
}

/* ------------------------------------------------------------------ */
/*  isortii: sort int64 array ix, optionally carry int64 array iy     */
/* ------------------------------------------------------------------ */
void isortii_(ITG *ix, ITG *iy, ITG *n, ITG *kflag)
{
    ITG nn = *n;
    if (nn < 2) return;

    int descending = (*kflag < 0);
    int carry      = (*kflag == 2 || *kflag == -2);

    lsd_sort_i64(ix, carry ? iy : NULL, NULL, nn);

    if (descending) {
        reverse_i(ix, nn);
        if (carry) reverse_i(iy, nn);
    }
}

/* ------------------------------------------------------------------ */
/*  isortid: sort int64 array ix, optionally carry double array dy    */
/* ------------------------------------------------------------------ */
void isortid_(ITG *ix, double *dy, ITG *n, ITG *kflag)
{
    ITG nn = *n;
    if (nn < 2) return;

    int descending = (*kflag < 0);
    int carry      = (*kflag == 2 || *kflag == -2);

    lsd_sort_i64(ix, NULL, carry ? dy : NULL, nn);

    if (descending) {
        reverse_i(ix, nn);
        if (carry) reverse_d(dy, nn);
    }
}

/* ------------------------------------------------------------------ */
/*  dsort: sort double array dx, optionally carry int64 array iy      */
/*                                                                     */
/*  Key-only: direct qsort on dx.                                     */
/*  With carry: pack into {double, ITG} structs, qsort, unpack.       */
/* ------------------------------------------------------------------ */

typedef struct { double key; ITG val; } DIPair;

static int cmp_double_asc(const void *a, const void *b)
{
    double da = *(const double *)a;
    double db = *(const double *)b;
    return (da > db) - (da < db);
}
static int cmp_double_desc(const void *a, const void *b)
{
    return cmp_double_asc(b, a);
}
static int cmp_dipair_asc(const void *a, const void *b)
{
    double da = ((const DIPair *)a)->key;
    double db = ((const DIPair *)b)->key;
    return (da > db) - (da < db);
}
static int cmp_dipair_desc(const void *a, const void *b)
{
    return cmp_dipair_asc(b, a);
}

void dsort_(double *dx, ITG *iy, ITG *n, ITG *kflag)
{
    ITG nn = *n;
    if (nn < 1) {
        fprintf(stderr,
                "*ERROR in dsort: the number of values to be sorted is "
                "not positive: %lld\n", (long long)nn);
        exit(201);
    }

    ITG kk = (*kflag < 0) ? -(*kflag) : *kflag;
    if (kk != 1 && kk != 2) {
        fprintf(stderr,
                "*ERROR in dsort: the sort control parameter is "
                "not 2, 1, -1 or -2\n");
        exit(201);
    }

    if (nn < 2) return;

    int descending = (*kflag < 0);
    int carry      = (*kflag == 2 || *kflag == -2);

    if (!carry) {
        qsort(dx, (size_t)nn, sizeof(double),
              descending ? cmp_double_desc : cmp_double_asc);
        return;
    }

    /* Pack, sort, unpack */
    ITG i;
    DIPair *pairs = (DIPair *)malloc((size_t)nn * sizeof(DIPair));
    for (i = 0; i < nn; i++) { pairs[i].key = dx[i]; pairs[i].val = iy[i]; }
    qsort(pairs, (size_t)nn, sizeof(DIPair),
          descending ? cmp_dipair_desc : cmp_dipair_asc);
    for (i = 0; i < nn; i++) { dx[i] = pairs[i].key; iy[i] = pairs[i].val; }
    free(pairs);
}
