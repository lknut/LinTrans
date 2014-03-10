namespace lawa {

template <typename T>
Option1D<T,Put>::Option1D()
    : optionparameters(),wr(0.0),useLinearTransformation(false), dwdM(-1)

{

}

template <typename T>
Option1D<T,Put>::Option1D(const OptionParameters1D<T,Put> &_optionparameters)
    : optionparameters(_optionparameters),wr(0.0),useLinearTransformation(false),dwdM(-1)

{
    assert(optionparameters.strike>=0.);
    assert(optionparameters.maturity>=0.);
    singularPoints.engine().resize(1);
    singularPoints(1) = 0.;//100;//0.;
	//LK this->calc_w();
}

template <typename  T>
Option1D<T,Put>::Option1D( const Option1D<T,Put>& _option )
	: optionparameters(_option.optionparameters), singularPoints(_option.singularPoints),values(_option.values),wr(_option.wr),dwdG(_option.dwdG),useLinearTransformation(_option.useLinearTransformation),dwdM(_option.dwdM),dwdr(_option.dwdr)
{
};


	//LK swap function that enables the use of a copy constructor 	
template <typename T>
void
Option1D<T,Put>::swap(Option1D<T,Put>& rhs){
		std::swap(optionparameters,rhs.optionparameters);
		std::swap(singularPoints,rhs.singularPoints);
		std::swap(values,rhs.values);
		std::swap(wr,rhs.wr);
		std::swap(dwdG,rhs.dwdG);
		std::swap(dwdM,rhs.dwdM);
		std::swap(dwdr,rhs.dwdr);
		std::swap(useLinearTransformation,rhs.useLinearTransformation);
}
	
template <typename T>	
Option1D<T,Put> Option1D<T,Put>::operator= (Option1D<T,Put> const& rhs)
{
  if (this != &rhs)  //oder if (*this != rhs)
  {
    /* kopiere elementweise, oder:*/
	Option1D<T,Put> tmp(rhs); //Copy-Konstruktor
	this->swap(tmp);
	//return tmp;
  }
  return *this; //Referenz auf das Objekt selbst zurückgeben
}

template <typename T>
T
Option1D<T,Put>::payoff(T S) const
{
	if (useLinearTransformation){
		return std::max( optionparameters.strike - std::exp( wr*optionparameters.maturity+std::log(S)), (T)0.);
	} else return std::max( optionparameters.strike - S, (T)0.);
}

template <typename T>
T
Option1D<T,Put>::payoff_log(T x) const
{
	if (useLinearTransformation){
		return optionparameters.strike*std::max(1. - std::exp( wr*optionparameters.maturity+x), (T)0.);
	} else return optionparameters.strike*std::max(1.-std::exp(x), (T)0.);
		
    
}


//payoff derivative with respect to model parameters
template <typename T>
T
Option1D<T,Put>::payoff_log(T x, ProcessType1D _processtype) const
{
	//strike=S_0!!!
	//std::cout << " dwdG = " << dwdG << std::endl;
	assert(useLinearTransformation);
	
	T mat = optionparameters.maturity;
	T strike = optionparameters.strike;
	T derivative_w=0.0;
	    // Sensitivities G
	if (_processtype==CGMYeSensitivitiesG){
		derivative_w = dwdG;
		//Sensitivities M
	}else if (_processtype==CGMYeSensitivitiesM){
		derivative_w = dwdM;
		//Sensitivities r
	}else if (_processtype==CGMYeSensitivitiesr){
		derivative_w = 1.0;
	} else {
		return payoff_log(x);
	}
	
	if (x > - wr*mat){
		return 0;
	}
	else {
		return  - derivative_w * mat * strike * std::exp( wr*mat + x  ) ;
	}
	
	
			 	
    
}

	




template <typename T>
T
Option1D<T,Put>::value(const ProcessParameters1D<T,BlackScholes> &processparameters, T S, T t) const
{
    assert(optionparameters.earlyExercise==false);
    static boost::math::normal norm;
    T r        = processparameters.r;
    T sigma    = processparameters.sigma;
    T strike   = optionparameters.strike;
    T maturity = optionparameters.maturity;
    T d_1 = ( log(S/strike) + ( r + sigma*sigma/2.0 ) * (maturity-t) )
                 / ( sigma * sqrt(maturity-t) );
    T d_2 = d_1 - sigma * sqrt(maturity-t);
    return strike*exp(-r*(maturity-t))*cdf(norm, -d_2) - S * cdf(norm, -d_1);
}

template <typename T>
T
Option1D<T,Put>::value(const ProcessParameters1D<T,CGMY> &processparameters, T S, T t)
{
    typedef typename std::map<T,T>::const_iterator const_it;

    assert(optionparameters.earlyExercise==false);
    assert(processparameters.k_Y!=0);
    assert(processparameters.k_Y!=1);
    CharacteristicFunction1D<T,CGMY> charfunc(processparameters);

    FourierPricer1D<T, CGMY> fp(charfunc, S, optionparameters.maturity,
                                std::max(optionparameters.strike-10.,(T)0.),
                                optionparameters.strike+10.);

    if (t==0) {
        const_it it_option = values.find(S);
        if (it_option!=values.end()) return (*it_option).second;
        else {
            fp.solve(10000,20);
            T val = fp(optionparameters.strike);
            values[S] = val;
            return val;
        }
    }
    else {
        //fp.solve(10000,17);
        fp.solve(10000,20);
        //fp.solve(10000,23);
    }
    return fp(optionparameters.strike);
}

template <typename T>
T
Option1D<T,Put>::value(const ProcessParameters1D<T,CGMYe> &processparameters, T S, T t)
{
	//std::cout << " in value " << std::endl;
	//bool tmpLinTrans = useLinearTransformation;
	//useLinearTransformation=false;
    typedef typename std::map<T,T>::const_iterator const_it;

    assert(optionparameters.earlyExercise==false);
    assert(processparameters.k_Y!=0);
    assert(processparameters.k_Y!=1);
    CharacteristicFunction1D<T,CGMYe> charfunc(processparameters);

    FourierPricer1D<T, CGMYe> fp(charfunc, S, optionparameters.maturity-t,
                                 std::max(optionparameters.strike-10.,(T)0.1),
                                 optionparameters.strike+10.);

    if (t==0) {
        const_it it_option = values.find(S);
        if (it_option!=values.end()) return (*it_option).second;
        else {
            std::cerr << "   -> Option1D<T,Put>: No high precision is used for computation of reference"
                      << " values!" << std::endl;
            //fp.solve(10000,16);
            fp.solve(10000,20);
            T val = fp(optionparameters.strike);
            values[S] = val;
            return val;
        }
    }
    else {
        //fp.solve(40000,17);
        fp.solve(10000,20);
        //fp.solve(10000,20);
        //fp.solve(10000,16);
    }
	//useLinearTransformation=tmpLinTrans;
	//std::cout << " useLinearTransformation = " << useLinearTransformation << std::endl;
    return fp(optionparameters.strike);
}

//LK, here we calculate the w of the Tops2012 Model
template <typename T>
void
Option1D<T,Put>::calc_w( const  ProcessParameters1D<T,CGMYe>  &_processparameters)
{
	useLinearTransformation=true;
		
	T C = _processparameters.k_C;
	T G = _processparameters.k_G;
	T M = _processparameters.k_M;
	T Y = _processparameters.k_Y;
						
	//calculating w from MScThesis_Tops, Appendix B
	// calculating w_1^1=w_i^j
	T w1;
	w1 = C*boost::math::tgamma(-Y)*( std::pow(M- 1,Y)-std::pow(M,Y)+Y*std::pow(M,Y-1)*1
									 +  std::pow(G+ 1,Y)-std::pow(G,Y)-Y*std::pow(G,Y-1) );
	wr = -w1;

	dwdG =- C*boost::math::tgamma(-Y)*(Y* std::pow(G+ 1,Y-1)- Y*std::pow(G,Y-1)-Y*(Y-1)*std::pow(G,Y-2) );

	dwdM =- C*boost::math::tgamma(-Y)*( Y*std::pow(M- 1,Y-1)- Y*std::pow(M,Y-1)+Y*(Y-1)*std::pow(M,Y-2));

	dwdr = 1;
	
	//LK Sigma influences the drift part as well
	T s = _processparameters.sigma*_processparameters.sigma;
	T drift   =  -0.5* s ; 
		
	//LK adding diffusion drift
	wr += drift;
	//LK adding risk-free interest rate
	wr += _processparameters.r;
	//LK adding singularPoint of the payoff function
	singularPoints(1) = - (wr*optionparameters.maturity);
}

	
}   // namespace lawa
