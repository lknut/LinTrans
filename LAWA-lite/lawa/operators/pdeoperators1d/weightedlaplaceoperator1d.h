/*
  LAWA - Library for Adaptive Wavelet Applications.
  Copyright (C) 2008-2013  Sebastian Kestler, Mario Rometsch, Kristina Steih, 
  Alexander Stippler.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 */


#ifndef LAWA_OPERATORS_PDEOPERATORS1D_WEIGHTEDLAPLACEOPERATOR1D_H
#define LAWA_OPERATORS_PDEOPERATORS1D_WEIGHTEDLAPLACEOPERATOR1D_H 1

#include <lawa/settings/enum.h>
#include <lawa/functiontypes/function.h>
#include <lawa/methods/adaptive/datastructures/index.h>
#include <lawa/integrals/integral.h>

namespace lawa {

/* Weighted Laplace OPERATOR 1D
 *
 *    a(v,u) =  Integral(w(x) * v_x * u_x)
 *
 */
template <typename T, typename Basis, QuadratureType Quad>
class WeightedLaplaceOperator1D{
    
    public:
        
        const Basis& basis;

        WeightedLaplaceOperator1D(const Basis& _basis, Function<T> weightFct);

        T
        operator()(XType xtype1, int j1, long k1,
                   XType xtype2, int j2, long k2) const;

        T
        operator()(const Index1D &row_index, const Index1D &col_index) const;
        
    private:
    
        Function<T> W;

        IntegralF<Quad, Basis, Basis> integral;

};

} // namespace lawa

#include <lawa/operators/pdeoperators1d/weightedlaplaceoperator1d.tcc>

#endif  // LAWA_OPERATORS_PDEOPERATORS1D_WEIGHTEDLAPLACEOPERATOR1D_H

