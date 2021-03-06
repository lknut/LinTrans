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

#ifndef APPLICATIONS_FINANCE_FOURIERPRICING_CGMYECHARACTERISTICFUNCTION1D_H
#define APPLICATIONS_FINANCE_FOURIERPRICING_CGMYECHARACTERISTICFUNCTION1D_H 1

#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <boost/math/special_functions/gamma.hpp>

#include <applications/finance/fourierpricing/characteristicfunction1d.h>
#include <applications/finance/processes/processes.h>


namespace lawa {

template <typename T>
struct CharacteristicFunction1D<T,CGMYe>
{
    CharacteristicFunction1D(const ProcessParameters1D<T,CGMYe> &_processparameters);

	//LK, neu -> todo Beschreibung
	CharacteristicFunction1D(T _r,T _C1, T _C2, T _C, T _G, T _M, T _Y, T _sigma1, T _sigma2, T _sigma, T _drift);
	
	
    gsl_complex
    operator()(T t, T v);

    gsl_complex
    operator()(T t, gsl_complex z);

    gsl_complex
    zeta_0();

    ProcessParameters1D<T,CGMYe> processparameters;
    T drift;                    //asserts that the underlying discounted process is martingale!
};

}   // namespace lawa

#include <applications/finance/fourierpricing/cgmyecharacteristicfunction1d.tcc>

#endif  // APPLICATIONS_FINANCE_FOURIERPRICING_CGMYECHARACTERISTICFUNCTION1D_H
