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

#ifndef APPLICATIONS_FINANCE_OPERATORS_SENSITIVITIES_CGMYEOPERATOR1D_R_H
#define APPLICATIONS_FINANCE_OPERATORS_SENSITIVITIES_CGMYEOPERATOR1D_R_H 1

#include <lawa/constructions/constructions.h>
#include <lawa/operators/deltas.h>
#include <lawa/settings/settings.h>
//#include <applications/finance/kernels/cgmykernel_G.h>
#include <applications/finance/operators/financeoperator1d.h>
#include <applications/finance/processes/processes.h>
//LK, new
#include <applications/finance/kernels/sensitivities/cgmykernel_G.h>

namespace lawa {

template <typename T, typename Basis1D>
struct FinanceOperator1D<T, CGMYeSensitivitiesr, Basis1D>
{
    #ifdef TRONE
        typedef typename std::tr1::unordered_map<Entry<Index1D>, T, entry_hashfunction<Index1D>,
                                                                 entry_eqfunction<Index1D> > EntryMap;
    #elif BOOST
        typedef typename boost::unordered_map<Entry<Index1D>, T, entry_hashfunction<Index1D>,
                                                                     entry_eqfunction<Index1D> > EntryMap;
    #else
        typedef typename __gnu_cxx::hash_map<Entry<Index1D>, T, entry_hashfunction<Index1D>,
                                                                 entry_eqfunction<Index1D> > EntryMap;
    #endif
    typedef typename EntryMap::value_type val_type;


    FinanceOperator1D(const Basis1D& _basis, const ProcessParameters1D<T,CGMYe> &_processparameters,
                      T _R1=0., T _R2=1., int order=10,
                      int internal_compression_level=-1
                      //T _convection=0., T _reaction=0.,
                      //const bool _use_predef_convection=true, const bool _use_predef_reaction=true
					  );

    void
    setCompressionLevel(int _internal_compression_level);

    T
    operator()(XType xtype1, int j1, int k1, XType xtype2, int j2, int k2) const;

    T
    operator()(const Index1D &row_index, const Index1D &col_index) const;

    const Basis1D                                    &basis;
    //const ProcessParameters1D<T,CGMYe>               &processparameters;
    //Kernel<T,CGMYeSensitivitiesG>                     kernel;
	//LK tmp
	//Kernel<T,CMGY>                                   kernel_tmp;
	//T                                                R1, R2;
    int                                              internal_compression_level;
    //T                                                convection, reaction;
    //const bool                                       use_predef_convection, use_predef_reaction;
    //T                                                OneDivSqrtR2pR1, OneDivR2pR1, R1DivR1pR2;

    //Integral<Gauss, Basis1D, Basis1D>                integral;
	//LK, has to be changed to CGMYeSensitivitiesG at some point of time!!!
	//SingularIntegral<Kernel<T,CGMY>,Basis1D,Basis1D> singularIntegral;
    //SingularIntegral<Kernel<T,CGMYeSensitivitiesG>,Basis1D,Basis1D> singularIntegral;

    mutable std::map<T,T> values_tailintegral;

    mutable EntryMap data;

};

}   // namespace lawa

#include <applications/finance/operators/sensitivities/cgmyeoperator1d_r.tcc>

#endif  // APPLICATIONS_FINANCE_OPERATORS_CGMYEOPERATOR1D_R_H
