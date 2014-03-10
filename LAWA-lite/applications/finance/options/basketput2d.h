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


#ifndef APPLICATIONS_FINANCE_OPTIONS_BASKETPUT2D_H
#define APPLICATIONS_FINANCE_OPTIONS_BASKETPUT2D_H 1

#include <boost/math/distributions/normal.hpp>

#include <applications/finance/options/option.h>
#include <applications/finance/options/optionparameters2d.h>
#include <applications/finance/processes/processes.h>
#include <applications/finance/fourierpricing/fourierpricing.h>

#include <applications/finance/mcpricing/nd_generator2d.hpp>
#include <applications/finance/mcpricing/cgmy.hpp>


#include <boost/math/special_functions/gamma.hpp> //LK: to include gamma for w-calculation!
#include <vector> // LK, for saving cgmy simulations

//test
//#include <iostream>

namespace lawa {

template <typename T>
struct Option2D<T,BasketPut>
{
    Option2D(void);

    Option2D(const OptionParameters2D<T,BasketPut> &_optionparameters);

	T
    payoff(T S1, T S2, T _w1, T _w2) const;

    T
    payoff(T S1, T S2) const;

    T
    payoff_log(T x1, T x2) const;

    T
    value(const ProcessParameters2D<T,BlackScholes2D> &processparameters, T S1, T S2, T t);

	//LK: processparameters include Sigma, Matrix of linear transformation, if processparam.lineartransformation=true
    T
    value(const ProcessParameters2D<T,CGMYeUnivariateJump2D> &processparameters, T S1, T S2, T t);

   T
   value(const ProcessParameters2D<T,CGMYeLinearTransformation2DJump2D> &processparameters, T S1, T S2, T t);

	
    void
    setNumberOfMCRuns(int _N);

	void
	setLevyMeasureTruncationEps(double _eps);	

	T
	get_w1();

	T
	get_w2();

	void
	calc_w();

	OptionParameters2D<T,BasketPut> optionparameters;
	
private:
    std::map<std::pair<T,T>,T> values;
    int N;
	double eps;
	double Tmt; //LK maturity - t, to be able to reuse simulations
	T w1,w2; //LK, weight from Tops-Model
	//  the following data structures save cgmy approximations
	vector<T> cgmyPlusDiffusion1,cgmyPlusDiffusion2; 
	//bool alte_methode; //LK tmp
	double eps_old;
	double dwda1;
	double dwda2;

};

} // namespace lawa

#include <applications/finance/options/basketput2d.tcc>

#endif // APPLICATIONS_FINANCE_OPTIONS_BASKETPUT2D_H
