namespace lawa {

 template <typename T>
 Option2D<T,BasketPut>::Option2D(void)
     : optionparameters(), N(300000),eps(0.000001),Tmt(-1.),eps_old(-1.),w1(-1.0)

 {
	 
 }

template <typename T>
Option2D<T,BasketPut>::Option2D(const OptionParameters2D<T,BasketPut> &_optionparameters)
    : optionparameters(_optionparameters), N(300000),eps(0.000001),Tmt(-1.),eps_old(-1.)

{
	// assertion is that we use the CGMYe model!
    assert(optionparameters.strike>=0.);
    assert(optionparameters.maturity>=0.);
	//if (optionparameters.lineartransformation){
	this->calc_w();
}


	//LK, TODO
	//
	// calculating w+r in one function and adjusting Sigma
	// then we don't need to store the processparameters
	//
	//
	// and Copy-Constructor and =-Operator implementation!
	
	

template <typename T>
T
Option2D<T,BasketPut>::payoff(T S1, T S2, T _w1, T _w2) const
{
	// OUTPUT: payoff(e^{ (r+w)(T) +_w + \Sigma Y_{T-t} }) = payoff(e^{ (r+w)(T-t) + \Sigma Y_{T-t} })
	// so _w is supposed to be -(r+w)*t
	//std::cout << " payoff called, lineartransformation= " << optionparameters.lineartransformation << std::endl;
	if (optionparameters.lineartransformation){
		T r= optionparameters.processparameters.r;
		T maturity = optionparameters.maturity;
		
		T SigmaS1= std::exp( (r + w1) * maturity
							 + _w2
							 + optionparameters.processparameters.Sigma11 * std::log(S1)
							 + optionparameters.processparameters.Sigma12 * std::log(S2));
		T SigmaS2= std::exp( (r + w2) * maturity
							 + _w2
							 + optionparameters.processparameters.Sigma21 * std::log(S1)
							 + optionparameters.processparameters.Sigma22 * std::log(S2));
				 
		return std::max(optionparameters.strike - optionparameters.weight1*SigmaS1 - optionparameters.weight2*SigmaS2, (T)0.);
	}
	else {
		return std::max(optionparameters.strike - optionparameters.weight1*S1 - optionparameters.weight2*S2, (T)0.);
	}

}


template <typename T>
T
Option2D<T,BasketPut>::payoff(T S1, T S2) const
{
	return payoff(S1,S2,0.0,0.0);
}

template <typename T>
T
Option2D<T,BasketPut>::payoff_log(T x1, T x2) const
{
	//LK, this payoff is only right at maturity -> no momentary payoff possible
	//LK, only holds if S_0 = strike
	
	// Sigma = [ 1- alpha, alpha ; alpha, 1- alpha]
	if (optionparameters.lineartransformation){
		T r= optionparameters.processparameters.r;
		T maturity = optionparameters.maturity;
		T Sigmax1= (r + w1)* maturity + optionparameters.processparameters.Sigma11 * x1 + optionparameters.processparameters.Sigma12 * x2;
		T Sigmax2= (r + w2)* maturity + optionparameters.processparameters.Sigma21 * x1 + optionparameters.processparameters.Sigma22 * x2;
		// std::cout << "LinearTransformation in basketput2d.tcc"
		// 		  << "Sigmax1 = " << Sigmax1 << "Sigmax2 = " << Sigmax2
		// 		  <<std::endl;

		// actual payoff, not for sensitivities
		// LK 
		double payoff = optionparameters.strike * std::max(1.-optionparameters.weight1*std::exp(Sigmax1)
												-optionparameters.weight2*std::exp(Sigmax2), (T)0.);
		// This is needed for the initial condition for sensitivities with respect to alpha
		double payoff1 = -optionparameters.weight1 * optionparameters.strike * ( maturity * dwda1 + (x2-x1)) * std::exp( Sigmax1  );
		double payoff2 = -optionparameters.weight2 * optionparameters.strike * ( maturity * dwda2 + (x1-x2)) * std::exp( Sigmax2  ); 

		return payoff;
		// This is the payoff for the initial condition for sensitivities with respect to alpha
		//if ( payoff > 0 ) return payoff1 + payoff2;
		//else return 0;
	}
	else {
		return optionparameters.strike*std::max(1.-optionparameters.weight1*std::exp(x1)
												-optionparameters.weight2*std::exp(x2), (T)0.);
	}
}


template <typename T>
T
Option2D<T,BasketPut>::value(const ProcessParameters2D<T,BlackScholes2D> &processparameters,
                             T S1, T S2, T t)
{
    typedef typename std::map<std::pair<T,T>,T>::const_iterator const_it;
    std::pair<T,T> pairS(S1,S2);

    const_it it = values.find(pairS);


	
    if (it != values.end()) {
        std::cerr << "Values found." << std::endl;
        return (*it).second;
    }
    else {
        T q11      = processparameters.sigma1*processparameters.sigma1;
        T q22      = processparameters.sigma2*processparameters.sigma2;
        T q12      = processparameters.sigma1*processparameters.sigma2*processparameters.rho;
        T q21      = q12;
        T drift1   = -0.5*q11; //LK: only if maturity=1?
        T drift2   = -0.5*q22;
        T r        = processparameters.r;
        T maturity = optionparameters.maturity;
        ND_Generator2D<T>               nd_generator2d(q11,q12,q21,q22,drift1,drift2);

        T optionPrice_estim = 0.;
        for (int j=1; j<=N; ++j) {
            T x1 = 0., x2 = 0.;
            nd_generator2d(maturity,x1,x2);
            T val = payoff(S1*exp(r*maturity+x1),S2*exp(r*maturity+x2));
            optionPrice_estim += exp(-r*maturity)*val;
        }
        optionPrice_estim /= N;
        values[pairS] = optionPrice_estim;
        std::cout.precision(16);
        std::cout << "Option price: " << optionPrice_estim << std::endl;
        return optionPrice_estim;
    }
}
	
	
	

template <typename T>
T
Option2D<T,BasketPut>::value(const ProcessParameters2D<T,CGMYeLinearTransformation2DJump2D> &processparameters,T S1, T S2, T t)
{
	// returns V(t,S1,S2)
	// S1,S2 are supposed to be
	// which equals e^{-r (T-t)} E[e^{(r+w)(T-t)+\Sigma Y_{T-t} }  ]
	typedef typename std::map<std::pair<T,T>,T>::const_iterator const_it;
	std::pair<T,T> pairS(S1,S2);

	const_it it = values.find(pairS);

	if (it != values.end()) {
        std::cerr << "Values found." << std::endl;
        return (*it).second;
    }

	//LK, this is for MC-Simulation -> pretend that Linear Transformation is used, but restore settings at the end
	//ProcessParameters2D<T,CGMYeUnivariateJump2D> tmp_processparamters(optionparameters.processparameters);
	bool tmpLinTrans = !optionparameters.lineartransformation;
	if (tmpLinTrans){
		optionparameters.processparameters=processparameters;
		optionparameters.lineartransformation=true;
		this->calc_w();
	}

	
	T maturity = optionparameters.maturity;
	//LK cgmy-Simulation reuse possible, but simulated paths will be the same (of course) -> dependencies (not wanted)
	bool simulation_reuse=true;
	if ((cgmyPlusDiffusion1.size() < N) || (Tmt != (maturity-t)) || (eps_old!=eps) ){
		cgmyPlusDiffusion1.resize(N);
		cgmyPlusDiffusion2.resize(N);
		simulation_reuse=false;
		Tmt = maturity -t;
		eps_old = eps;
	}
	else{
		simulation_reuse=true;
	}


	// read from file ------------------------------------------------------------------------------
	// this is used to save cgmy + gaussian simulations in order to obtain a fast (quasi-optimal) MC
	// calculation
	bool save_sim_to_file=false;
	std::stringstream simulation_file;
	simulation_file << "mcpricing/saved_simulations/Sim_T1_N" << N  <<"_eps" << eps <<"_paths" << processparameters.getSimulationName() << ".csv";

	std::ofstream simfile;

	//LK, reading already simulated and saved values
	std::ifstream read_simfile;
	read_simfile.open(simulation_file.str().c_str());
	 	 
	if( ( (!read_simfile.good())  ) && (!simulation_reuse)){
		simfile.open(simulation_file.str().c_str());
		//LK, the simulated values will be saved
		save_sim_to_file=true;
		std::cout << " Saving paths of the processes " << std::endl;
	} else if(!simulation_reuse){
		cgmyPlusDiffusion1.resize(N);
		std::string str;
		int i;
		for (i=1; ((i<=N) && (!read_simfile.eof())); ++i){
			if (i % 100000 == 0) std::cout << "MC, read simulations, i = " << i << " of " << N << std::endl;
			getline(read_simfile,str,' ');
			cgmyPlusDiffusion1[i-1] = atof(str.c_str());
			getline(read_simfile,str,' ');
			cgmyPlusDiffusion2[i-1]=atof(str.c_str());
			simulation_reuse=true; //LK: use vector data
		}
		if (i<N){
			std::cout << "Sim-File not correct (not enough lines) -> Abort, i = " << i << ", N = " << N << std::endl;
			exit(2);
		}
		// file has not enough simulations (maybe it is still written into)
	}	 
	// end reading from file ----------------------------------------------------------------------


	T C1 = processparameters.k_C1, C2 = processparameters.k_C2;
	T G1 = processparameters.k_G1, G2 = processparameters.k_G2;
	T M1 = processparameters.k_M1, M2 = processparameters.k_M2;
	T Y1 = processparameters.k_Y1, Y2 = processparameters.k_Y2;
	T Sigma11 = processparameters.Sigma11, Sigma12 = processparameters.Sigma12;
	T Sigma21 = processparameters.Sigma21, Sigma22 = processparameters.Sigma22;
	T r        = processparameters.r;
	
	//LK simulate independent CGMY processes with center zero!
	// test if e^{-rT} e{(r+w)T + X_t} is a martingale!!
	T exp_estimator1 = 0.0;
	T exp_estimator2 = 0.0;

	//double 
	CGMYSimulator csim1(C1,G1,M1,Y1,eps);
	CGMYSimulator csim2(C2,G2,M2,Y2,eps);

	T cgmy1,cgmy2;
		
	T q11      = processparameters.sigma1*processparameters.sigma1;
	T q22      = processparameters.sigma2*processparameters.sigma2;
	T q12      = processparameters.sigma1*processparameters.sigma2*processparameters.rho;
	T q21      = q12;

	
	ND_Generator2D<T> nd_generator2d(q11,q12,q21,q22,(T)0,(T)0); //drift1,drift2);
	T x1,x2;


	T val;
	T optionPrice_estim = 0.;

	for (int j=1; j<=N; ++j) {
		if ((j % 500==0) && (simulation_reuse==false)) std::cout <<" MC Simulation "<< std::setw(6) << j  << "/" << N
																 << ", MC-Approx of V(" << t <<","<<S1<<","<<S2<<") = "
																 << optionPrice_estim/(j-1)
																 <<  std::endl;
		if (!simulation_reuse){
			x1 = 0.;
			x2 = 0.;
			nd_generator2d(maturity-t,x1,x2);
			cgmy1 = csim1.simZeroCenter(maturity-t);
			cgmy2 = csim2.simZeroCenter(maturity-t);
			//CGMY and Diffusion saved together
			cgmyPlusDiffusion1[j-1]= cgmy1 + x1;
			cgmyPlusDiffusion2[j-1]= cgmy2 + x2;
			if (save_sim_to_file){
				simfile << cgmy1 + x1 << " " << cgmy2 + x2<<" " <<std::endl;
			}
			
		} else{
			x1=0.;
			x2=0.;
			// CGMY and Diffusion saved together
			cgmy1 = cgmyPlusDiffusion1[j-1];
			cgmy2 = cgmyPlusDiffusion2[j-1];
		}
			
		//exp_estimator1 += std::exp(w1*maturity + x1 +cgmy1);
		//exp_estimator2 += std::exp(w2*maturity + x2 +cgmy2);

		// funktioniert nur, da \Sigma = \Id!!
		
		// -(r+w1) funktioniert nur bei Sigma=Id!!, S1,S2 sind untransformierte
		// neue payoff funktion? payoff(S1,S2 (to transform), w1,w2 (let untransfomred)
		// und payoff(S1,S2) = payoff(S1,S2,0,0)?
		val = payoff(  S1*exp(-+x1+cgmy1)  ,  S2*exp( +x2+cgmy2)  , - (r + w1)*t,  - (r + w2)*t  );
		optionPrice_estim += exp(-r*(maturity-t))*val;
	}
	if (!simulation_reuse){
		simfile.close();
	}	
	optionPrice_estim /= N;
	values[pairS] = optionPrice_estim;
	
	//restore settings if LinTrans-setting was just for MC-Sim
	if(tmpLinTrans){
		//LK they are not used anyway
		//optionparameters.processparameters=tmp_processparameters;
		optionparameters.lineartransformation=false;
		w1=0.0;
		w2=0.0;
	}
	return optionPrice_estim;
}


template <typename T>
T
Option2D<T,BasketPut>::value(const ProcessParameters2D<T,CGMYeUnivariateJump2D> &_processparameters,T S1, T S2, T t)
{
	T r = _processparameters.r;
	T sigma1 = _processparameters.sigma1;
	T sigma2 = _processparameters.sigma2;
	T rho = 0.;

	T C1 = _processparameters.k_C1, C2= _processparameters.k_C2;
	T G1 = _processparameters.k_G1, G2= _processparameters.k_G2;
	T M1 = _processparameters.k_M1, M2= _processparameters.k_M2;
	T Y1 = _processparameters.k_Y1, Y2= _processparameters.k_Y2;

	
	//LK: Transformation Matrix, max(||) <! G,M
	T Sigma11=1.0, Sigma12=0.0, Sigma21=0.0, Sigma22=1.0;
	
	ProcessParameters2D<T,CGMYeLinearTransformation2DJump2D> linTransProcessparameters(r, sigma1, sigma2, rho,
																					   C1,  G1, M1, Y1,
																					   C2,  G2, M2, Y2,
																					   Sigma11, Sigma12,
																					   Sigma21, Sigma22);

	return this->value(linTransProcessparameters,S1,S2,t, true);
}
	

template <typename T>
void
Option2D<T,BasketPut>::setNumberOfMCRuns(int _N) {
    N = _N;
}

template <typename T>
void
Option2D<T,BasketPut>::setLevyMeasureTruncationEps(double _eps) {
    eps = _eps;
}



template <typename T>
T
Option2D<T,BasketPut>::get_w1(){
		return w1;
};

template <typename T>
T
Option2D<T,BasketPut>::get_w2(){
		return w2;
};


template <typename T>
void
Option2D<T,BasketPut>::calc_w()
{
	if (optionparameters.lineartransformation==true){
		//LK: we need an independent model for the LinTransform case, but this will be stored generally
		assert(optionparameters.processparameters.rho==0.);
		//LK, parameters needed	
		T C1 = optionparameters.processparameters.k_C1, C2 = optionparameters.processparameters.k_C2;
		T G1 = optionparameters.processparameters.k_G1, G2 = optionparameters.processparameters.k_G2;
		T M1 = optionparameters.processparameters.k_M1, M2 = optionparameters.processparameters.k_M2;
		T Y1 = optionparameters.processparameters.k_Y1, Y2 = optionparameters.processparameters.k_Y2;
			
		T Sigma11 = optionparameters.processparameters.Sigma11, Sigma12 = optionparameters.processparameters.Sigma12;
		T Sigma21 = optionparameters.processparameters.Sigma21, Sigma22 = optionparameters.processparameters.Sigma22;
	
						
		//calculating w from MScThesis_Tops, Appendix B
		// calculating w_1^1=w_i^j
		T w11,w12;
		w11 = C1*boost::math::tgamma(-Y1)*( std::pow(M1- Sigma11,Y1)-std::pow(M1,Y1)+Y1*std::pow(M1,Y1-1)*Sigma11
											+  std::pow(G1+ Sigma11,Y1)-std::pow(G1,Y1)-Y1*std::pow(G1,Y1-1)*Sigma11 );
		w12 = C2*boost::math::tgamma(-Y2)*( std::pow(M2- Sigma12,Y2)-std::pow(M2,Y2)+Y2*std::pow(M2,Y2-1)*Sigma12
											+  std::pow(G2+ Sigma12,Y2)-std::pow(G2,Y2)-Y2*std::pow(G2,Y2-1)*Sigma12 );
		w1 = -w11 - w12;
		
		//LK calculating w_2
		T w21,w22;
		w21 = C1*boost::math::tgamma(-Y1)*( std::pow(M1- Sigma21,Y1)-std::pow(M1,Y1)+Y1*std::pow(M1,Y1-1)*Sigma21
											+  std::pow(G1+ Sigma21,Y1)-std::pow(G1,Y1)-Y1*std::pow(G1,Y1-1)*Sigma21 );
		w22 = C2*boost::math::tgamma(-Y2)*( std::pow(M2- Sigma22,Y2)-std::pow(M2,Y2)+Y2*std::pow(M2,Y2-1)*Sigma22
											+  std::pow(G2+ Sigma22,Y2)-std::pow(G2,Y2)-Y2*std::pow(G2,Y2-1)*Sigma22 );
		w2 = -w21 - w22;

		//LK Sigma influences the drift part as well
		T s1 = optionparameters.processparameters.sigma1*optionparameters.processparameters.sigma1; //sigma_1²
		T s2 = optionparameters.processparameters.sigma2*optionparameters.processparameters.sigma2; // sigma_2²
		T drift1   = -0.5* ( Sigma11 * Sigma11 * s1 + Sigma12 * Sigma12 * s2   ); 
		T drift2   = -0.5* (  Sigma21 * Sigma21 * s1 + Sigma22 * Sigma22 * s2  );
		
		//LK adding diffusion drift
		w1 += drift1;
		w2 += drift2;



		//LK, this is needed to calculate the sensitivities with respect to alpha, where
		// \Sigma = [ (1-alpha) alpha; alpha (1- alpha) ]
		// note: this is temporarily done, alpha should not be fixed here!
		double alpha = Sigma12;
		// LK, we have to make sure that Sigma is symmetric, and that Sigma11=Sigma22
		// ji=11
		double dwda11 = C1* boost::math::tgamma(-Y1) * ( 1* Y1* std::pow((M1 - (1-alpha)  ),Y1-1)  - Y1*std::pow(M1,Y1-1)
														 - Y1* std::pow(G1+ (1-alpha),Y1-1) + 	Y1*std::pow(G1,Y1-1) );
			//ji = 12	
		double dwda12 = C2* boost::math::tgamma(-Y2) * ( - Y2* std::pow((M2 - alpha ),Y2-1)  + Y2*std::pow(M2,Y2-1)
														 + Y2* std::pow(G2+ alpha,Y2-1) - 	Y2*std::pow(G2,Y2-1) );
			//ji = 11
		double dbda11 = - s1 * (1- alpha) ;
		double dbda12 =  s2 * alpha;
		dwda1 = - ( dwda11 + dbda11 + dwda12 + dbda12);


		// ji = 21
		double dwda21 = C1* boost::math::tgamma(-Y1) * ( - Y1* std::pow(M1 - alpha,Y1-1)  + Y1*std::pow(M1,Y1-1)
														 + Y1* std::pow(G1+ alpha,Y1-1) -  	Y1*std::pow(G1,Y1-1) );
			//ji = 22	
		double dwda22 = C2* boost::math::tgamma(-Y2) * ( Y2*std::pow(M2 - (1-alpha),Y2-1)  - Y2*std::pow(M2,Y2-1)
														 - Y2* std::pow(G2+ (1-alpha),Y2-1) + 	Y2*std::pow(G2,Y2-1) );
			//ji = 21
		double dbda21 =  s1 * alpha ;
		double dbda22 =  s2 * (-1) * (1- alpha);
		dwda2 = - ( dwda21 + dbda21 + dwda22 + dbda22);

		std::cout << " we calculate sensitivities with respect to alpha!! " << std::endl;
		std::cout << " dwda1 = "  << dwda1 << std::endl;
		std::cout << " dwda2 = " << dwda2 << std::endl;
				
				
	}
}

}   //namespace lawa
