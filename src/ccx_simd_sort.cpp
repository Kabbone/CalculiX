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
 * SIMD-accelerated replacements for selected Fortran sort subroutines.
 *
 * Replaces: isortii.f  isortid.f  dsort.f
 *
 * Uses the x86-simd-sort header-only library (Intel, Apache-2.0).
 * Runtime ISA dispatch: SSE/AVX2/AVX-512 selected automatically.
 *
 * Fortran calling convention (gfortran):
 *   - symbol name = lowercase subroutine name + trailing underscore
 *   - all arguments passed by pointer
 *   - integers are int64_t throughout (-fdefault-integer-8)
 *
 * kflag semantics (same as original SLATEC routines):
 *    1  sort key ascending,  ignore payload
 *   -1  sort key descending, ignore payload
 *    2  sort key ascending,  carry payload
 *   -2  sort key descending, carry payload
 */

#include "x86simdsort.h"

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <vector>

/* Fortran integer type: always 8 bytes with -fdefault-integer-8 */
typedef long long ITG;

extern "C" {

/* ------------------------------------------------------------------ */
/* isortii: sort int64 key array ix, optionally carry int64 array iy  */
/*                                                                     */
/* ix(*)  - key array                                                  */
/* iy(*)  - payload array (used only for kflag = ±2)                  */
/* n      - number of elements                                         */
/* kflag  - control flag (see header)                                  */
/* ------------------------------------------------------------------ */
void isortii_(ITG *ix, ITG *iy, ITG *n, ITG *kflag)
{
    ITG nn = *n;
    if (nn < 1) return;

    bool descending = (*kflag < 0);
    bool carry      = (*kflag == 2 || *kflag == -2);

    if (carry) {
        x86simdsort::keyvalue_qsort(ix, iy, nn, false, descending);
    } else {
        x86simdsort::qsort(ix, nn, false, descending);
    }
}

/* ------------------------------------------------------------------ */
/* isortid: sort int64 key array ix, optionally carry double array dy  */
/*                                                                     */
/* ix(*)  - key array                                                  */
/* dy(*)  - double payload array (used only for kflag = ±2)           */
/* n      - number of elements                                         */
/* kflag  - control flag (see header)                                  */
/* ------------------------------------------------------------------ */
void isortid_(ITG *ix, double *dy, ITG *n, ITG *kflag)
{
    ITG nn = *n;
    if (nn < 1) return;

    bool descending = (*kflag < 0);
    bool carry      = (*kflag == 2 || *kflag == -2);

    if (!carry) {
        x86simdsort::qsort(ix, nn, false, descending);
        return;
    }

    /*
     * x86-simd-sort's keyvalue_qsort requires key and value to have the
     * same width (both 4-byte or both 8-byte).  int64+double are both
     * 8 bytes so the template instantiation is valid.
     *
     * If a future version of the library changes this constraint, fall
     * back to the argsort + gather approach:
     *
     *   auto idx = x86simdsort::argsort(ix, nn, false, descending);
     *   std::vector<ITG>    tmp_k(nn);
     *   std::vector<double> tmp_v(nn);
     *   for (ITG i = 0; i < nn; i++) { tmp_k[i]=ix[idx[i]]; tmp_v[i]=dy[idx[i]]; }
     *   std::copy(tmp_k.begin(), tmp_k.end(), ix);
     *   std::copy(tmp_v.begin(), tmp_v.end(), dy);
     */
    x86simdsort::keyvalue_qsort(ix, dy, nn, false, descending);
}

/* ------------------------------------------------------------------ */
/* dsort: sort double key array dx, optionally carry int64 array iy   */
/*                                                                     */
/* dx(*)  - double key array                                           */
/* iy(*)  - int64 payload array (used only for kflag = ±2)            */
/* n      - number of elements                                         */
/* kflag  - control flag (see header)                                  */
/* ------------------------------------------------------------------ */
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

    bool descending = (*kflag < 0);
    bool carry      = (*kflag == 2 || *kflag == -2);

    if (!carry) {
        x86simdsort::qsort(dx, nn, false, descending);
        return;
    }

    /*
     * double key + int64 value: same width (8+8 bytes), valid template.
     * See the fallback comment in isortid_ if this ever needs replacing.
     */
    x86simdsort::keyvalue_qsort(dx, iy, nn, false, descending);
}

} /* extern "C" */
