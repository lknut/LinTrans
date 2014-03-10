#include <fstream>
#include <iostream>
#include <flens/flens.h>
#include <lawa/lawa.h>
#include <applications/finance/fourierpricing/fourierpricing.h>
#include <applications/finance/initialconditions/initialconditions.h>
#include <applications/finance/options/options.h>
#include <applications/finance/operators/operators.h>
#include <applications/finance/processes/processes.h>
#include <applications/finance/righthandsides/righthandsides.h>
//LK, neu
#include <applications/finance/operators/sensitivities/cgmyeoperator1d_G.h>
#include <applications/finance/operators/sensitivities/cgmyeoperator1d_M.h>
#include <applications/finance/operators/sensitivities/cgmyeoperator1d_r.h>
//LK, test
#include <lawa/methods/uniform/algorithms/assembler1d.h>
#include <lawa/methods/uniform/solvers/thetascheme1d_LTI_Sensitivities.h>
//#include <lawa/methods/uniform/solvers/thetascheme1d_LTI_Sensitivities.h>

//#define SOLVER_DEBUG 1

using namespace std;
using namespace lawa;

///  Typedefs for precision. Important: `long double` is only available for $L_2$-
///  orthonormal constructions!
//LK, was typedef /*long*/ double T;
typedef /*long*/ double T;

///  Typedefs for Flens data types:
typedef flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >  FullColMatrixT;
typedef flens::SparseGeMatrix<flens::CRS<T,flens::CRS_General> >    SparseMatrixT;
typedef flens::DiagonalMatrix<T>                                    DiagonalMatrixT;
typedef flens::DenseVector<flens::Array<T> >                        DenseVectorT;

///  Set up of the option parameters
const OptionType1D optiontype = Put;
T strike = 1.;
T maturity = 1.;
T S0 = 1.;


//For error measuremeant
T delta;

//LK to store output folder
string outputFolder;

///  Definition of the process type
const ProcessType1D  processtype  = CGMYe;
//LK, new
const ProcessType1D  sensitivitytype = CGMYeSensitivitiesG;
//const ProcessType1D  sensitivitytype = CGMYeSensitivitiesG;
//const ProcessType1D  sensitivitytype = CGMYeSensitivitiesr;
//const ProcessType1D  processtype  = BlackScholes;
//const ProcessType1D  processtype  = CGMY;

/* Reference values for Europ. option from Fang & Oosterlee (2008) and Almendral & Oosterlee (2007)
 * Option parameters strike = 100, maturity = 1, SO = 100
 * Process parameters for CGMY: r = 0.1, C = 1, G = M = 5, Y =1.5   49.790905469 (call)
 * Process parameters for CGMY: r = 0.1, C = 1, G = M = 5, Y =0.8   14.789424 (put)
 * Process parameters for CGMY: r = 0.1, C = 1, G = M = 5, Y =0.1    6.353404 (put)
 */

///  Typedefs for basis construction:
//typedef Basis<T,Primal,Interval,Dijkema>                      Basis1D;
typedef Basis<T,Orthogonal,Interval,Multi>                      Basis1D;

///  Integral definitions:
typedef Integral<Gauss,Basis1D,Basis1D>                         IntegralBasis1DBasis1D;
typedef IntegralF<Gauss,Basis1D>                                IntegralFBasis1D;

///  Operator definitions. The first one is required for the implementation of the $\theta$-scheme,
///  the second one is the realization of the CGMY (or, as a special case, the Black-Scholes)
///  operator.
typedef IdentityOperator1D<T, Basis1D>                          ScalarProduct1D;
typedef FinanceOperator1D<T, processtype, Basis1D>              FinanceOp;
//LK new
typedef FinanceOperator1D<T, sensitivitytype, Basis1D>          SensitivitiesOp;

///  Definition for the right-hand side required for the barrier option approach.
typedef OptionRHS1D<T, optiontype, processtype, Basis1D>        OptionRhs;

///  Definition of the time stepping method
typedef ThetaScheme1D_LTI<T, Basis1D, FinanceOp, OptionRhs>     ThetaStepScalarProduct1D;
// theta-stepping for sensitivities 
typedef ThetaScheme1D_LTI_Sensitivities<T, Basis1D,
										FinanceOp,
										SensitivitiesOp,
										OptionRhs>              ThetaStepSensitivities1D;

typedef TimeStepping<T, ThetaStepScalarProduct1D>               TimeStepperScalarProduct1D;

///  Routine for the computation of the option pricing error arising from localization and
///  discretization in space and time.
template<typename T, OptionType1D OType, ProcessType1D PType, typename Basis>
void
ComputeL2ErrorAndPlotSolution(Option1D<T,OType> &option,
                              ProcessParameters1D<T,PType> &processparameters,
                              const Basis &basis, int J, const DenseVectorT &u,const DenseVectorT &tilde_u, const DenseVectorT &tilde_u0,
							  const DenseVectorT &tilde_u_compr,
                              const DenseVectorT &u0, T R1, T R2, bool excessToPayoff,
                              T &L2error, T &Linftyerror, T &Linftyerror_sensitivities);

///  Computation of the wavelet basis coefficients of the approximation of the initial condition
///  (only required for barrier option approach).
void
getPu0(const Basis1D &basis, DenseVectorT &Pu0,  const Option1D<T,Put> &option, T R1, T R2, int J, ProcessType1D _processtype,  ProcessParameters1D<T,CGMYe> &processparameters);

template<typename T, OptionType1D OType, ProcessType1D PType>
void
plotOptionPriceSurface(Option1D<T,OType> &option, const ProcessParameters1D<T,PType> &processparameters, T R1, T R2);

template<typename T, ProcessType1D PType, typename Basis>
void
spyStiffnessMatrix(const Basis &basis, T R1, T R2, const FinanceOp &finance_op, int J, T tol,
                   ProcessParameters1D<T,PType> &processparameters, bool pattern);

int
main(int argc, char *argv[])
{
    cout.precision(12);
    if ((argc < 14) || (argc>15) ) { // " 1  2 3  4  5     6         7 8     9 101112    13       14
		                                     //R1 R2 excessToPayoff theta timesteps r sigma G M Y delta sencmpr? outputFolder(opt)" << endl;
        cerr << "usage: " << argv[0] << " j0 J R1 R2 theta timesteps r sigma G M Y delta sencmpr? outputFolder(opt)" << endl;
        exit(1);
    }
    ///  Wavelet basis parameters
    int d=2, d_=2;
	//LK: sensitivities only work for d=2 so far!
    int j0    = atoi(argv[1]);
    int j_max = atoi(argv[2]);
	delta = T(atof(argv[12]));
	outputFolder= "";
	if (argc==15){
		outputFolder=argv[14];
		if (outputFolder[outputFolder.length()-1] != '/') outputFolder= outputFolder +'/';//strcat(outputFolder,"/");
	}

    ///  Localization parameters: Considered interval is $[-R_1,R_2]$ and the localization approach
    ///  `sencmpr? = 0` (...) and `sencmor? = 1` (...).
    T   R1         = T(atof(argv[3]));
    T   R2         = T(atof(argv[4]));
    //int etp        = atoi(argv[5]);
	int sens        = atoi(argv[13]);

    ///  Parameters for the $\theta$-time stepping method.
    T theta        = T(atof(argv[5]));
    int timesteps  = atoi(argv[6]);

    ///  Parameters for the CGMY financial model
    T   r     = T(atof(argv[7]));
    T   sigma = T(atof(argv[8]));
    T   G     = T(atof(argv[9]));
    T   M     = T(atof(argv[10]));
    T   Y     = T(atof(argv[11]));

    //bool excessToPayoff   = (etp == 1) ? true : false;
	bool sensitivitiesCmpr = (sens == 1) ? true : false;
    //ProcessParameters1D<T,BlackScholes>   processparameters(0.04, 0.2);

    //ProcessParameters1D<T,CGMY>           processparameters(r, 1., G, M, Y);
    //ProcessParameters1D<T,CGMY>           processparameters(0.1, 1., 5., 5., 1.5 );
    //ProcessParameters1D<T,CGMY>           processparameters(0.1, 1., 5., 5., 0.8 );
    //ProcessParameters1D<T,CGMY>           processparameters(0.1, 1., 5., 5., 0.1 );
    //ProcessParameters1D<T,CGMY>           processparameters(0.04, 1., 2.4, 4.5, 1.8 );

    ///  Initialization of the CGMY model (here, extended with an additional diffusion term)
    ProcessParameters1D<T,CGMYe>          processparameters(r, 1., G, M, Y, sigma);
    //ProcessParameters1D<T,CGMYe>          processparameters(0.04, 1., 2.4, 4.5, 1.8, 0.1 );
    //ProcessParameters1D<T,CGMYe>          processparameters(0.04, 1., 7.4, 8.5, 1.1, 0.1 );
    //ProcessParameters1D<T,CGMYe>          processparameters(0.04, 1., 5., 5., 1.5, 0.1 );
    cout << processparameters << endl;

	//OptionParameters1D<T,Put> optionparameters(strike, maturity, false);
	OptionParameters1D<T,Put> optionparameters(strike, maturity, false);  //, true, processparameters);

    ///  Note: We exclude explicit schemes.
    if (theta < 0.5) {
        cout << "theta should be larger than 0.5!" << endl;
        exit(1);
    }

    ///  Number of employed time steps
    T                       timestep = optionparameters.maturity/T(timesteps);

    ///  Basis initialization
    //LK, was Basis1D   basis(d,d_,j0); // for Dijkema basis
    Basis1D   basis(d,j0);  // for $L_2$-orth. basis
    basis.enforceBoundaryCondition<DirichletBC>();

    int                             order=20;
    Option1D<T,optiontype>  option(optionparameters);
	
	//LK: this is needed to use the Linear transformation model -> different payoff function
	option.calc_w(processparameters);

	//plotOptionPriceSurface(option, processparameters, R1, R2);
    //return 0;

    std::stringstream filename;
    /*if (excessToPayoff) {
        filename << outputFolder << "wavelet_galerkin_option_pricer1d_conv_" << d << "_" << d_
                 << "_R1_" << R1 << "_R2_" << R2 << "_wetp__"
                 << strike << "_" << maturity << "_" << processparameters << ".txt";
				 } */
    //else {
    filename << outputFolder << "wavelet_galerkin_option_pricer1d_conv_" << d << "_" << d_
             << "_R1_" << R1 << "_R2_" << R2 << "_woetp__"
             << strike << "_" << maturity << "_" << processparameters << ".txt";
		//}
    std::ofstream convfile(filename.str().c_str());
	
    for (int J=j0; J<=j_max; ++J) {

        Timer time;
        time.start();

        ///  Wavelet basis coefficients for the solution and the initial condition.
        DenseVectorT u(basis.mra.rangeI(J)), u0(basis.mra.rangeI(J)), tilde_u(basis.mra.rangeI(J)),
			// two vectors (u_compr,tilde_u_compr) needed for a for compressed rhs of the sensitivity PIDE approach of HilberSchwabWinter2008   
			tilde_u0(basis.mra.rangeI(J)),u_compr(basis.mra.rangeI(J)), tilde_u_compr(basis.mra.rangeI(J))  ;

		
        /// approximation of the initial condition
        getPu0(basis, u0, option, R1, R2, J,processtype, processparameters);
        // approximation of the sensitivity intiial condition!
		getPu0(basis, tilde_u,option,R1,R2,J,sensitivitytype, processparameters);

		// initial condition is the same for a compressed rhs
		tilde_u0=tilde_u;
		tilde_u_compr=tilde_u;
				
        ///  Initialization of the operator associated to the financial model
		//LK                                       int_compr_level, _conv,_reac,_use_prefed_convection=true, _use_predef_reaction=true)
        FinanceOp                               finance_op(basis, processparameters, R1, R2, order,
														   J,0.0,0.0,!option.useLinearTransformation);
		//LK Sensitivity Operator Initialization (with compression)
		SensitivitiesOp                         sensitivities_opCompr(basis, processparameters, R1, R2,
																	  order, J);//, 0.0,0.0,false); // J was -1.0
		//LK Sensitivity Operator Initialization (without compression)
	    SensitivitiesOp                         sensitivities_op(basis, processparameters, R1, R2,
																 order, -1.0);//, 0.0,0.0,false); // J was -1.0


        ///  Initialization of the right-hand side vector for a general method.
        OptionRhs                         rhs(optionparameters, processparameters, basis,
                                              R1, R2, 0); // was excessToPayoff);

        ///  Initialization of the solver for each time step
        ThetaStepScalarProduct1D          scheme(theta, basis, finance_op, rhs,
                                                 true, false, 0., 1e-12);
		ThetaStepSensitivities1D          scheme_sensitivities(theta,basis, finance_op, sensitivities_op,rhs,
															   true, false, 0., 1e-12);
		ThetaStepSensitivities1D          scheme_sensitivitiesCompr(theta,basis, finance_op, sensitivities_opCompr,rhs,
															   true, false, 0., 1e-12);

		
        ///  Initialization of the solver for the whole time stepping
        TimeStepperScalarProduct1D        timestepmethod(scheme, timestep, timesteps, J);

        ///  Solving the fully discretized system
        //LK, was
		//u = timestepmethod.solve(u0, false);
		 u=u0;
		 //tmp!!!
		 //u=tilde_u0;
		 u_compr=u0;
		   //std::cout << "tilde_u = " << std::endl << tilde_u << std::endl;
		 for(int k = 1; k <= timesteps; ++k){
			 if (k%50==0) std::cout << "level " << J << ", time step = " << k << " of " << timesteps << std::endl;
			 if (sensitivitiesCmpr){
				 //std::cout << " ok" << std::endl;
				 scheme_sensitivities.solve((k-1)*timestep, k*timestep, u, tilde_u, J); //levelX);
			 }
			 scheme_sensitivitiesCompr.solve((k-1)*timestep, k*timestep, u_compr, tilde_u_compr, J); //levelX);
			 // 	//std::cout << "k = " << k << ": u_next = " << u_next << std::endl;
		 }
		
        time.stop();
        T comp_time = time.elapsed();
		std::cout << " duration: " << comp_time << " sec" << std::endl;

        T L2error = 0.;
        T Linftyerror = 0.;
		T Linftyerror_sensitivities=0.;
		//LK, abdiskont neu std::exp(-r * maturity)*
		T approx = std::exp(-r * maturity)*(1./sqrt(R1+R2))*evaluate(basis, J, u_compr, (0+R1)/(R1+R2), 0);
		
        ComputeL2ErrorAndPlotSolution(option, processparameters, basis, J, u_compr, tilde_u, tilde_u_compr, tilde_u0, u0, R1, R2,
                                      false, L2error, Linftyerror, Linftyerror_sensitivities);
		//           1                   2                        3 
        cout      << J << " " << basis.mra.cardI(J) << " " << comp_time << " "
		//             4                   5           
				  << L2error << " " << Linftyerror << " value(S_0) = "<< approx << endl;
		//           1                     2                     3
        convfile  << J << " " << basis.mra.cardI(J) << " " << comp_time << " "
			//         4                  5                            6	
                  << L2error << " " << Linftyerror <<  " " << Linftyerror_sensitivities << endl;
		//spyStiffnessMatrix(basis, R1, R2, finance_op,J , 0.0000000000001,processparameters, false);
    }

 

    return 0;
}

template<typename T, OptionType1D OType, ProcessType1D PType, typename Basis>
void
ComputeL2ErrorAndPlotSolution(Option1D<T,OType> &option,
                              ProcessParameters1D<T,PType> &processparameters,
                              const Basis &basis, int J,
							  const DenseVectorT &u,
							  const DenseVectorT &tilde_u,
							  const DenseVectorT &tilde_u_compr,
							  const DenseVectorT &tilde_u0,
                              const DenseVectorT &u0, T R1, T R2, bool excessToPayoff,
                              T &L2error, T &Linftyerror, T& Linftyerror_sensitivities)
{
    std::stringstream filename;
    filename << outputFolder << "solutionPlot_" << J <<".txt";
    std::ofstream plotFile(filename.str().c_str());


    TruncatedPutOption1D<T,OType> truncatedPutOption;
    truncatedPutOption.setOption(option); //,processparameters);
	
    T h1 = (R1+R2)*pow2i<T>(-7);
    truncatedPutOption.setTruncation(-R1,R2,1,h1);

    T tmp_maturity = option.optionparameters.maturity;
    T tmp_strike   = option.optionparameters.strike;
    T r        = processparameters.r;

    /// LK: Skalierung? 1/(sqrt...) da die wavelets sonst nicht L2-bdd sind? sqrt, da L2 = \int (\cdot)^2 ~ 1
    T approxPutS0;
	if (option.useLinearTransformation){
		approxPutS0= exp(-r*tmp_maturity)*(1./sqrt(R1+R2))*evaluate(basis, J, u, (R1)/(R1+R2), 0);
	} else approxPutS0= exp(-r*tmp_maturity)*(1./sqrt(R1+R2))*evaluate(basis, J, u, (r*tmp_maturity+R1)/(R1+R2), 0);
    T approxCallS0 = approxPutS0 + S0 - tmp_strike*exp(-processparameters.r*tmp_maturity);
	T exactPutS0   = option.value(processparameters, tmp_strike, 0);
    T exactCallS0  = exactPutS0 + S0 - tmp_strike*exp(-processparameters.r*tmp_maturity);

    std::cout << "Reference value put : " << exactPutS0 << ", approximated value: " << approxPutS0 << std::endl;
    std::cout << "Reference value call: " << exactCallS0 << ", approximated value: " << approxCallS0 << std::endl;

    L2error     = 0.;
    Linftyerror = 0.;
	// measures L_{\infty}-error of the sensitivites of the uncompressed tildeA-scheme vs. the compressed version
	Linftyerror_sensitivities=0.;
    //LK, was T h = pow2i<T>(-J-2)*(R1+R2);
	//T h = pow2i<T>(-J-2)*(R1+R2);
	T h = pow2i<T>(-8-4)*(R1+R2);
    //T h=0.03;
    //T h = pow2i<T>(-J-2)*(R1+R2);
    //T delta = 1.;
    std::cerr << "Error computation started..." << std::endl;
    for (T x=-delta*R1; x<=delta*R2; x+=h) {
		//LK: from LinTrans.cpp
		//T approx = std::exp(-r*maturity)* (1./sqrt(R1+R2))*evaluate(basis, J, u, (x-r*maturity+R1)/(R1+R2), 0);
		//
        //LK, was todo $std::exp((-r) * maturity)*$
		//std::cout << " r = " << r << " maturity = " << maturity << std::endl;
		T disc = std::exp((-r) * maturity);
		//tmp, without disc* !
		T approx = disc* (1./sqrt(R1+R2))*evaluate(basis, J, u, (x+R1)/(R1+R2), 0);
		//LK, check if discount necessary!!
        T approx_u0 = (1./sqrt(R1+R2))*evaluate(basis, J, u0, (x+R1)/(R1+R2), 0);
		T approx_tilde_u =disc* (1./sqrt(R1+R2))*evaluate(basis, J, tilde_u, (x+R1)/(R1+R2), 0);
		T approx_tilde_u_compr = disc* (1./sqrt(R1+R2))*evaluate(basis, J, tilde_u_compr, (x+R1)/(R1+R2), 0);
		if (sensitivitytype == CGMYeSensitivitiesr){
			approx_tilde_u+= - maturity*  approx;
			approx_tilde_u_compr+= - maturity *  approx;   
		}
		T approx_tilde_u0 = (1./sqrt(R1+R2))*evaluate(basis, J, tilde_u0, (x+R1)/(R1+R2), 0);
        T exact = 0.;
		T spot;
        //LK, was: T spot = tmp_strike*std::exp(x-r*tmp_maturity);
		if (option.useLinearTransformation){
			spot= tmp_strike*std::exp(x);
		} else spot = tmp_strike*std::exp(x-r*tmp_maturity);
        exact = option.value(processparameters, spot, 0);
		//exact = option.value(processparameters, spot, 0);

        if (excessToPayoff) exact -= option.payoff(tmp_strike*exp(x));
        //if (excessToPayoff) approx += option.payoff(tmp_strike*exp(x));

        if ((fabs(x+delta*R1)<1e-12) || (fabs(x-delta*R2) < 1e-12)) {
            L2error += 0.5*h*std::pow(approx-exact,(T)2.);
        }
        else    {
            L2error += h*std::pow(approx-exact,(T)2.);
        }
        Linftyerror = std::max(Linftyerror, fabs(exact-approx));
		
		Linftyerror_sensitivities = std::max(Linftyerror_sensitivities, fabs(approx_tilde_u-approx_tilde_u_compr));
		//           1           2               3               4                      5
        plotFile  << x << " " << spot << " "<<  exact  << " " << approx <<  " " << option.payoff_log(x-option.wr*tmp_maturity) << " "
			//              6                                  7                         8
                  << truncatedPutOption.g_trunc(x) << " " << approx_u0   << " " <<  approx_tilde_u << " "
			//                                                 
			//	  << option.log_payoff_G_trunc(x)  <<" "
			//               9                                                    10
				  << truncatedPutOption.g_trunc(x,sensitivitytype)  <<" "  << approx_tilde_u0
			//                11                        
				  << " " << approx_tilde_u_compr << " " 
				  << endl;

    }
	//std::cout << "Linftyerror_Sensitivities = " << Linftyerror_sensitivities << std::endl;
    plotFile.close();
    L2error = sqrt(L2error);
    std::cerr << "...finished." << std::endl;
    //std::cerr << "Please hit enter..." << std::endl;
    //getchar();
   
}

void
getPu0(const Basis1D &basis, DenseVectorT &Pu0,  const Option1D<T,Put> &option, T R1, T R2, int J, ProcessType1D _processtype, ProcessParameters1D<T,CGMYe> &processparameters)
{
	//if (_processtype==processtype){
	//LK: ERROR!!!!!!!!!!!!!!!!!!!!!!! -> not right approximation
	// seems that wavelet-integrals are not right. scaling functions seem ok

	
	//LK, do I need truncation? -> initial condition has to be H1
	// => need this for non-zero initial-condition!
	TruncatedPutOption1D<T,optiontype> truncatedPutOption;
	truncatedPutOption.setOption(option);
	//LK was://T h = (R1+R2)*pow2i<T>(-7);
	T h = (R1+R2)*pow2i<T>(-8);
	//T h=0.;
	//LK, neu
	truncatedPutOption.setProcesstype(_processtype);
	//std::cout << " processtype = " << _processtype;
	truncatedPutOption.setTruncation(-R1,R2,1,h);

	//LK this is the initial function -> change to $\tilde{g}$
	//LK what is Function<T> ?
	Function<T> g_fct(truncatedPutOption.g_trunc_derivative,truncatedPutOption.singPts);
	
	InitialCondition1D<Basis1D> initialCondition1D(g_fct,basis,-R1,R2);
	int j0 = basis.j0;
	int N  = basis.mra.cardI(J);
	int pos = basis.mra.rangeI(J).firstIndex();
	for (int k=basis.mra.rangeI(j0).firstIndex(); k<=basis.mra.rangeI(j0).lastIndex(); ++k) {
		Pu0(pos) = initialCondition1D(XBSpline,j0,k,0);
		++pos;
	}
	for (int j=j0; j<J; ++j) {
		for (int k=basis.rangeJ(j).firstIndex(); k<=basis.rangeJ(j).lastIndex(); ++k) {
			Pu0(pos) = initialCondition1D(XWavelet,j,k,0);
			++pos;
		}
	}

	

}

template<typename T, OptionType1D OType, ProcessType1D PType>
void
plotOptionPriceSurface(Option1D<T,OType> &option, const ProcessParameters1D<T,PType> &processparameters, T R1, T R2)
{
    cout << "plotOptionPriceSurface called." << endl;

    std::stringstream filename;
    filename << "optionpricesurface_" << processparameters << ".txt";
    std::ofstream plotFile(filename.str().c_str());

    Kernel<T,CGMY> kernel(processparameters);

    T tmp_maturity = option.optionparameters.maturity;
    T tmp_strike   = option.optionparameters.strike;
    T r        = processparameters.r;

    T h = 0.2;
    for (T tau=0.1; tau<=tmp_maturity; tau+=0.1) {
        for (T x=-R1; x<=R2; x+=h) {
            plotFile << tau << " " << x << " "
                     << exp(r*tau)*option.value(processparameters,tmp_strike*exp(x-r*tau),tmp_maturity-tau) // kernel.ExpXmOnemX_k_pos+kernel.ExpXmOnemX_k_neg
                     << " " << option.payoff(tmp_strike*exp(x)) << endl;
        }
        plotFile << endl;
    }
}

template<typename T, ProcessType1D PType, typename Basis>
void
spyStiffnessMatrix(const Basis &basis, T R1, T R2, const FinanceOp &finance_op, int J, T tol,
                   ProcessParameters1D<T,PType> &processparameters, bool pattern)
{
    std::stringstream matrixfilename;
    matrixfilename << outputFolder << "A_" << J << "_" << basis.d << "_" << basis.d_ << "_" << tol
                   << "_R1_" << R1 << "_R2_" << R2 << processparameters << ".txt";
    int N =basis.mra.cardI(J);
    SparseMatrixT A(N,N);

    int j0 = basis.j0;
    int offsetJ = basis.rangeJ(j0).firstIndex()-1;
    int offsetI = basis.mra.rangeI(j0).firstIndex()-1;

    for(int k1 = basis.mra.rangeI(j0).firstIndex(); k1 <= basis.mra.rangeI(j0).lastIndex(); ++k1){
        for(int k2 = basis.mra.rangeI(j0).firstIndex(); k2 <= basis.mra.rangeI(j0).lastIndex(); ++k2){
            T val = finance_op(XBSpline, j0, k1, XBSpline, j0, k2);
            if(fabs(val) > tol){
                A(k1-offsetI, k2-offsetI) = pattern ? 1 : val;
            }
            else {
                cout << k1 << ", " << k2 << ": " << val << endl;
            }
        }
    }

    // SF * W
    for(int k1 = basis.mra.rangeI(j0).firstIndex(); k1 <= basis.mra.rangeI(j0).lastIndex(); ++k1){
        for(int j = j0; j <= J-1; ++j){
            for(int k2 = basis.rangeJ(j).firstIndex(); k2 <= basis.rangeJ(j).lastIndex(); ++k2){
                T val = finance_op(XBSpline, j0, k1, XWavelet, j, k2);
                if(fabs(val) > tol){
                    A(k1-offsetI,  basis.mra.cardI(j) + k2 - offsetJ) = pattern ? 1 : val;
                }
            }
        }
    }

    // W * SF
    for(int k2 = basis.mra.rangeI(j0).firstIndex(); k2 <= basis.mra.rangeI(j0).lastIndex(); ++k2){
        for(int j = j0; j <= J-1; ++j){
            for(int k1 = basis.rangeJ(j).firstIndex(); k1 <= basis.rangeJ(j).lastIndex(); ++k1){
                T val = finance_op(XWavelet, j, k1, XBSpline, j0, k2);
                if(fabs(val) > tol){
                    A(basis.mra.cardI(j) + k1 - offsetJ, k2 - offsetI) = pattern ? 1 : val;
                }
            }
        }
    }

        // W * W
    for(int j = j0; j <= J-1; ++j){
        for(int k1 = basis.rangeJ(j).firstIndex(); k1 <= basis.rangeJ(j).lastIndex(); ++k1){
            for(int j_ = j0; j_ <= J-1; ++j_){
                for(int k2 = basis.rangeJ(j_).firstIndex(); k2 <= basis.rangeJ(j_).lastIndex(); ++k2){
                    T val = finance_op(XWavelet, j, k1, XWavelet, j_, k2);
                    if(fabs(val) > tol){
                        A(basis.mra.cardI(j) + k1 - offsetJ, basis.mra.cardI(j_) + k2 - offsetJ) = pattern ? 1 : val;
                    }
                }
            }
        }
    }

    A.finalize();


    spy(A,matrixfilename.str().c_str(),true,(T)0.);
}


