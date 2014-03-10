/*
 *   Copyright (c) 2010, Michael Lehn
 *
 *   All rights reserved.
 *
 *   Redistribution and use in source and binary forms, with or without
 *   modification, are permitted provided that the following conditions
 *   are met:
 *
 *   1) Redistributions of source code must retain the above copyright
 *      notice, this list of conditions and the following disclaimer.
 *   2) Redistributions in binary form must reproduce the above copyright
 *      notice, this list of conditions and the following disclaimer in
 *      the documentation and/or other materials provided with the
 *      distribution.
 *   3) Neither the name of the FLENS development group nor the names of
 *      its contributors may be used to endorse or promote products derived
 *      from this software without specific prior written permission.
 *
 *   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 *   A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 *   OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 *   SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 *   LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 *   DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 *   THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 *   (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 *   OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef FLENS_BLAS_LEVEL3_R2K_H
#define FLENS_BLAS_LEVEL3_R2K_H

#include <cxxblas/cxxblas.h>
#include <flens/matrixtypes/general/impl/gematrix.h>
#include <flens/matrixtypes/hermitian/impl/hematrix.h>
#include <flens/matrixtypes/symmetric/impl/symatrix.h>

namespace flens { namespace blas {

//-- her2k
template <typename ALPHA, typename MA, typename MB, typename BETA, typename MC>
    void
    r2k(cxxblas::Transpose trans,
        const ALPHA &alpha,
        const GeMatrix<MA> &A, const GeMatrix<MB> &B,
        const BETA &beta,
        HeMatrix<MC> &C);

//-- syr2k
template <typename ALPHA, typename MA, typename MB, typename BETA, typename MC>
    void
    r2k(cxxblas::Transpose trans,
        const ALPHA &alpha,
        const GeMatrix<MA> &A, const GeMatrix<MB> &B,
        const BETA &beta,
        SyMatrix<MC> &C);

} } // namespace blas, flens

#include <flens/blas/level3/r2k.tcc>

#endif // FLENS_BLAS_LEVEL3_R2K_H