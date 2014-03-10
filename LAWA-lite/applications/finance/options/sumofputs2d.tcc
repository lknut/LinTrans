namespace lawa {

template <typename T>
Option2D<T,SumOfPuts>::Option2D(void)
    : optionparameters(), optionparameters1(), optionparameters2(),
      option1(), option2(), useLinearTransformation(false)

{

}

template <typename T>
Option2D<T,SumOfPuts>::Option2D(const OptionParameters2D<T,SumOfPuts> &_optionparameters)
    : optionparameters(_optionparameters),
      optionparameters1(_optionparameters.strike1, _optionparameters.maturity, _optionparameters.earlyExercise),
      optionparameters2(_optionparameters.strike2, _optionparameters.maturity, _optionparameters.earlyExercise),
      option1(optionparameters1), option2(optionparameters2),useLinearTransformation(false)

{
    assert(optionparameters.strike1>0.);
    assert(optionparameters.strike2>0.);
    assert(optionparameters.maturity>=0.);
}



	//LK, copy-Constructor
template <typename T>
Option2D<T,SumOfPuts>::Option2D( const Option2D<T,SumOfPuts>& _option )
	:   optionparameters(_option.optionparameters),
		optionparameters1(_option.optionparameters1),
		optionparameters2(_option.optionparameters2),
		option1(optionparameters1), option2(optionparameters2),
		useLinearTransformation(_option.useLinearTransformation),
		//singularPoints(_option.singularPoints), //LK, müsste klappen
		//values(_option.values), //Lk, müsste klappen
		wr1(_option.wr1),wr2(_option.wr2), //Lk, ok
		//option1D(_option.option1D), //Lk müsste klappen
		Sigma11(_option.Sigma11), Sigma12(_option.Sigma12),
		Sigma21(_option.Sigma21), Sigma22(_option.Sigma22){
};


	
//LK
template <typename T>
Option2D<T,SumOfPuts>
Option2D<T,SumOfPuts>::operator= (Option2D<T,SumOfPuts> const& rhs){
	if (this != &rhs)  //oder if (*this != rhs)
		{
			/* kopiere elementweise, oder:*/
			Option2D<T,SumOfPuts> tmp(rhs); //Copy-Konstruktor
			this->swap(tmp);
			//return tmp;
		}
	return *this; //Referenz auf das Objekt selbst zurückgeben
}



//LK
template <typename T>
void
Option2D<T,SumOfPuts>::swap(Option2D<T,SumOfPuts> & rhs){
	std::swap(optionparameters,rhs.optionparameters);
	std::swap(optionparameters1,rhs.optionparameters1);
	std::swap(optionparameters2,rhs.optionparameters2);
	//std::swap(singularPoints,rhs.singularPoints);
	//std::swap(values,rhs.values);
	std::swap(wr1,rhs.wr1);
	std::swap(wr2,rhs.wr2);
	//std::swap(dwdM,rhs.dwdM);
	//std::swap(dwdr,rhs.dwdr);
	std::swap(useLinearTransformation,rhs.useLinearTransformation);
	std::swap(Sigma11,rhs.Sigma11);
	std::swap(Sigma12,rhs.Sigma12);
	std::swap(Sigma21,rhs.Sigma21);
	std::swap(Sigma22,rhs.Sigma22);
};	


	
	//LK S1, S2 wrt independent coordinates Y	
template <typename T>
T
Option2D<T,SumOfPuts>::payoff(T S1, T S2, T tau) const
{

	T w1 = optionparameters.weight1;
	T w2 = optionparameters.weight2;
	T s1 = optionparameters.strike1;
	T s2 = optionparameters.strike2;
	T mat = optionparameters.maturity;
	T x1 = std::log(S1);
	T x2 = std::log(S2);

	if (useLinearTransformation){
		return w1*s1* std::max(1- std::exp( wr1 * mat + Sigma11 * x1 + Sigma12 * x2), (T)0.) +
			w2*s2* std::max(1- std::exp( wr2 * mat + Sigma21 * x1 + Sigma22 * x2), (T)0.);
	} else {
		std::cerr << " sumofputs2d.tcc Not using Linear Transformation!!! " << std::endl;
		return    w1* s1* std::max(1.-std::exp(x1),(T)0.)
           +  w2*s2*   std::max(1.-std::exp(x2),(T)0.);
	}
	
}		

	
template <typename T>
T
Option2D<T,SumOfPuts>::payoff(T S1, T S2) const
{
    return this->payoff(S1,S2,optionparameters.maturity);
		//optionparameters.weight1*std::max(optionparameters.strike1 - S1,(T)0.)
		//   + optionparameters.weight2*std::max(optionparameters.strike2 - S2,(T)0.);
}

	//LK, x independent coordinates	
template <typename T>
T
Option2D<T,SumOfPuts>::payoff_log(T x1, T x2) const
{
	T w1 = optionparameters.weight1;
	T w2 = optionparameters.weight2;
	T s1 = optionparameters.strike1;
	T s2 = optionparameters.strike2;
	T mat = optionparameters.maturity;
	//std::cout << "lt = " << useLinearTransformation << " Sigma11 = " << Sigma11 << "Sigma22 = " << Sigma22 << " w1 = " << w1 << " s1 = " << s1 << std::endl;
	if (useLinearTransformation){
		return w1*s1* std::max(1- std::exp( wr1 * mat + Sigma11 * x1 + Sigma12 * x2), (T)0.) +
			w2*s2* std::max(1- std::exp( wr2 * mat + Sigma21 * x1 + Sigma22 * x2), (T)0.);
	} else
		return    w1* s1* std::max(1.-std::exp(x1),(T)0.)
           +  w2*s2*   std::max(1.-std::exp(x2),(T)0.);
}

template <typename T>
T
Option2D<T,SumOfPuts>::value(const ProcessParameters2D<T,BlackScholes2D> &processparameters,
                             T S1, T S2, T t) const
{
    ProcessParameters1D<T,BlackScholes> processparameters1(processparameters.r, processparameters.sigma1);
    ProcessParameters1D<T,BlackScholes> processparameters2(processparameters.r, processparameters.sigma2);

    return    optionparameters.weight1*option1.value(processparameters1,S1,t)
            + optionparameters.weight2*option2.value(processparameters2,S2,t);
}

template <typename T>
T
Option2D<T,SumOfPuts>::value(const ProcessParameters2D<T,CGMYeUnivariateJump2D> &processparameters,
                             T S1, T S2, T t)
{
    return    optionparameters.weight1*option1.value(processparameters.proc_param1,S1,t)
            + optionparameters.weight2*option2.value(processparameters.proc_param2,S2,t);
}


	//////////////////////////////////
    // LK, Additional Test Value    //
    //////////////////////////////////


template <typename T>
T
Option2D<T,SumOfPuts>::value(const ProcessParameters2D<T,CGMYeLinearTransformation2DJump2D> &processparameters,
                             T _S1, T _S2, T t)
{

	assert(G1==G2);
	assert(M1==M2);
	assert(Y1 ==Y2);
	assert(Y1!=0);
	assert(Y1!=1);	
	T C1 = processparameters.k_C1, C2 = processparameters.k_C2;
	T G1 = processparameters.k_G1, G2 = processparameters.k_G2;
	T M1 = processparameters.k_M1, M2 = processparameters.k_M2;
	T Y1 = processparameters.k_Y1, Y2 = processparameters.k_Y2;
	T Sigma11 = processparameters.Sigma11, Sigma12 = processparameters.Sigma12;
	T Sigma21 = processparameters.Sigma21, Sigma22 = processparameters.Sigma22;
	T r        = processparameters.r;
	T maturity = optionparameters.maturity;


	T _C1 = (Sigma11 + Sigma21)*C1;
	T _C2 = (Sigma12 + Sigma22)*C2;
	//T C = _C1 + _C2;

	T q11      = processparameters.sigma1*processparameters.sigma1;
	T q22      = processparameters.sigma2*processparameters.sigma2;
	assert (_processparameters.rho==0);

	
	T sigma1 =std::sqrt( (Sigma11 *q11 + Sigma21 *q22));
	T sigma2 = std::sqrt( (Sigma21 *q11 + Sigma22 *q22) );
	//T sigma = std::sqrt((Sigma11 + Sigma21)* q11 +  (Sigma12 + Sigma22) * q22);

	ProcessParameters1D<T,CGMYe>  processparametersCGMYe1(r, _C1, G1,M1,Y1,sigma1);
	ProcessParameters1D<T,CGMYe>  processparametersCGMYe2(r, _C2, G2,M2,Y2,sigma2);

	T S1 = std::exp(Sigma11 * std::log(_S1) + Sigma12 * std::log(_S2));
	T S2 = std::exp(Sigma21 * std::log(_S1) + Sigma22 * std::log(_S2));
	
		
	return optionparameters.weight1 * option1.value(processparametersCGMYe1,S1,t) +
		optionparameters.weight2 * option2.value(processparametersCGMYe2,S2,t);
}
	

template <typename T>
void
Option2D<T,SumOfPuts>::calc_w(const ProcessParameters2D<T,CGMYeLinearTransformation2DJump2D> &processparameters){

	//std::cout <<" basketprodput2d.tcc, calc_w() " << std::endl;
	if (!useLinearTransformation){
		useLinearTransformation=true;
		T C1 = processparameters.k_C1, C2 = processparameters.k_C2;
		T G1 = processparameters.k_G1, G2 = processparameters.k_G2;
		T M1 = processparameters.k_M1, M2 = processparameters.k_M2;
		T Y1 = processparameters.k_Y1, Y2 = processparameters.k_Y2;
		assert(G1==G2);
		assert(M1==M2);
		assert(Y1 ==Y2);
		Sigma11 = processparameters.Sigma11;
		Sigma12 = processparameters.Sigma12;
		Sigma21 = processparameters.Sigma21;
		Sigma22 = processparameters.Sigma22;

		T r        = processparameters.r;
		T maturity = optionparameters.maturity;

		//calculating w from MScThesis_Tops, Appendix B
		// calculating w_1^1=w_i^j
		double w1, w11,w12;
		w11 = C1*boost::math::tgamma(-Y1)*( std::pow(M1- Sigma11,Y1)-std::pow(M1,Y1)+Y1*std::pow(M1,Y1-1)*Sigma11
											+  std::pow(G1+ Sigma11,Y1)-std::pow(G1,Y1)-Y1*std::pow(G1,Y1-1)*Sigma11 );
		w12 = C2*boost::math::tgamma(-Y2)*( std::pow(M2- Sigma12,Y2)-std::pow(M2,Y2)+Y2*std::pow(M2,Y2-1)*Sigma12
											+  std::pow(G2+ Sigma12,Y2)-std::pow(G2,Y2)-Y2*std::pow(G2,Y2-1)*Sigma12 );
		w1 = -w11 - w12;
		//LK calculating w_2
		double w2, w21,w22;
		w21 = C1*boost::math::tgamma(-Y1)*( std::pow(M1- Sigma21,Y1)-std::pow(M1,Y1)+Y1*std::pow(M1,Y1-1)*Sigma21
											+  std::pow(G1+ Sigma21,Y1)-std::pow(G1,Y1)-Y1*std::pow(G1,Y1-1)*Sigma21 );
		w22 = C2*boost::math::tgamma(-Y2)*( std::pow(M2- Sigma22,Y2)-std::pow(M2,Y2)+Y2*std::pow(M2,Y2-1)*Sigma22
											+  std::pow(G2+ Sigma22,Y2)-std::pow(G2,Y2)-Y2*std::pow(G2,Y2-1)*Sigma22 );
		w2 = -w21 - w22;
		
		std::cout << " w1 = " << w1 << " w2 = " << w2 << std::endl;

		//LK simulate independent CGMY processes with center zero!
		// double eps=0.000001;
		// CGMYSimulator csim1(processparameters.k_C1,processparameters.k_G1,processparameters.k_M1,
		// 					processparameters.k_Y1,eps);
		// CGMYSimulator csim2(processparameters.k_C2,processparameters.k_G2,processparameters.k_M2,
		// 					processparameters.k_Y2,eps);

		
		T q11      = processparameters.sigma1*processparameters.sigma1;
		T q22      = processparameters.sigma2*processparameters.sigma2;
		assert (processparameters.rho==0);
		T q12      = processparameters.sigma1*processparameters.sigma2*processparameters.rho;
		T q21      = q12;

		
		T drift1   = -0.5*(q11 * Sigma11*Sigma11 + Sigma12*Sigma12 * q22   ); 
		T drift2   = -0.5*(q11 * Sigma21*Sigma21 + Sigma22*Sigma22 * q22    );
		
		//LK adding diffusion drift
		w1 += drift1;
		w2 += drift2;

		std::cout <<" basketprodput2d.tcc, w1+drif = " << w1 << ", w2 + drift = " << w2 << std::endl;
		//LK, tmp!!!

		wr1 = w1+r;
		wr2 = w2 +r;

	}
	
}


template <typename T>
T
Option2D<T,SumOfPuts>::get_wr1(){
	return wr1; //- optionparameters.r;
}

template <typename T>
T
Option2D<T,SumOfPuts>::get_wr2(){
	return wr2; //- optionparameters.r;
}





	
	
}   //namespace lawa
