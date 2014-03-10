namespace lawa {


	/*
template <typename T, OptionType1D OType>
TruncatedPutOption1D<T, OType>::TruncatedPutOption()
	:processtype(CGMYe)
{

};*/	

template <typename T, OptionType1D OType>
ProcessType1D
TruncatedPutOption1D<T, OType>::processtype;
	
template <typename T, OptionType1D OType>
T
TruncatedPutOption1D<T, OType>::left;

template <typename T, OptionType1D OType>
T
TruncatedPutOption1D<T, OType>::right;

template <typename T, OptionType1D OType>
int
TruncatedPutOption1D<T, OType>::type;

template <typename T, OptionType1D OType>
T
TruncatedPutOption1D<T, OType>::truncWidth;

template <typename T, OptionType1D OType>
flens::DenseVector<Array<T> >
TruncatedPutOption1D<T, OType>::singPts;

template <typename T, OptionType1D OType>
Option1D<T,OType>
TruncatedPutOption1D<T, OType>::option;

template <typename T, OptionType1D OType>
void
TruncatedPutOption1D<T, OType>::setOption(const Option1D<T,Put> &_option)
{
	/* option.optionparameters = _option.optionparameters;
	option.values=_option.values;
	option.useLinearTransformation=_option.useLinearTransformation;
	option.wr=_option.wr;
	option.singularPoints   = _option.singularPoints;*/
	option=_option;
}


	
// 	//LK: for using the linear transformation model
// template <typename T, OptionType1D OType>
// void
// TruncatedPutOption1D<T, OType>::setOption(const Option1D<T,Put> &_option,
// 										  const  ProcessParameters1D<T,CGMYe>  &_processparamters)
// {
//     option.optionparameters = _option.optionparameters;
// 	//LK, new
// 	option.calc_w(_processparameters);
// 	option.singularPoints   = _option.singularPoints;
// }
	

template <typename T, OptionType1D OType>
void
TruncatedPutOption1D<T, OType>::setTruncation(T _left, T _right, int _type, T _truncWidth) 
{
	
    left       = _left;
    right      = _right;
    type       = _type;
    truncWidth = _truncWidth;
    int p = option.singularPoints.engine().length();
    singPts.engine().resize(p+1);
	singPts(1) = left+truncWidth;
    for (int i=1; i<=p; ++i) {
        singPts(1+i) = option.singularPoints(i);
	}
	
}

template <typename T, OptionType1D OType>
T
TruncatedPutOption1D<T, OType>::g_trunc(T x)
{
    if (x>=left+truncWidth) {
        return option.payoff_log(x);
    }
    else  {
        T b = left+truncWidth;
        T a = left;
        T val = option.payoff_log(b);
        return val*(x-a)/(b-a);
    }
}

	//LK: to obtain sensitivities

	//LK, new
template <typename T, OptionType1D OType>
T
TruncatedPutOption1D<T, OType>::g_trunc(T x, ProcessType1D _processtype)
{
    if (x>=left+truncWidth) {
		//LK schnell verbessert, sollte noch ordentlich angepasst werden
		//std::cout << " truncatedputoption1d.txx, ordentlich anpassen " << std::endl;
		if ( (processtype==CGMYeSensitivitiesG) || (processtype==CGMYeSensitivitiesM) || (processtype==CGMYeSensitivitiesr) ){
			if (std::abs(x-singPts(2)) <= truncWidth){
				T a = singPts(2)-truncWidth;
				T b = singPts(2) + truncWidth;
				T val = option.payoff_log(a, _processtype);
				return val * (b-x)/(b-a);
			}
		}
		return option.payoff_log(x, _processtype);
    }
    else  {
        T b = left+truncWidth;
        T a = left;
        T val = option.payoff_log(b, _processtype);
        return val*(x-a)/(b-a);
    }
}

	//LK, new
template <typename T, OptionType1D OType>
void
TruncatedPutOption1D<T, OType>::setProcesstype(ProcessType1D _processtype)
{
	processtype = _processtype;
}

	//LK, new
template <typename T, OptionType1D OType>
T
TruncatedPutOption1D<T, OType>::g_trunc_derivative(T x)
{
	// 8, 9, 10 sensitivity processes
	if  ( (processtype==CGMYeSensitivitiesG) || (processtype==CGMYeSensitivitiesM) || (processtype==CGMYeSensitivitiesr) ){
		//std::cout << " processtype = " << processtype << std::endl;
		//std::cout << " truncatedputoption1d.tcc if, processtype = " << processtype << std::endl;
		return g_trunc(x,processtype);
	} else {
		//std::cout  << std::endl<< " truncatedputoption1d.tcc, else processtype = " << processtype << std::endl;
		return g_trunc(x);
	}
	
	//return 0.0;
}	

}   // namespace lawa
