namespace lawa {

template <typename T>
OptionParameters1D<T,Put>::OptionParameters1D(void)
	: strike((T)0.), maturity((T)0.), earlyExercise(false)//, lineartransformation(false),processparameters()
{
	
}

template <typename T>
OptionParameters1D<T,Put>::OptionParameters1D(T _strike, T _maturity, bool _earlyExercise,bool _lineartransformation, ProcessParameters1D<T,CGMYe> _processparameters)
	: strike(_strike), maturity(_maturity), earlyExercise(_earlyExercise)//,
	  //lineartransformation(_lineartransformation),processparameters(_processparameters)
{
	//lineartransformation=true;
}

}   // namespace lawa
