#include <iostream>
#include <lawa/lawa.h>

#include <applications/finance/initialconditions/initialconditions.h>
#include <applications/finance/options/options.h>
#include <applications/finance/operators/operators.h>
#include <applications/finance/processes/processes.h>

//LK: string manipulation (outputFolder)
#include <string>


// definitions for programm use
//#define __CGMYeUnivariateJump2D
#define __CGMYeLinearTransformation2DJump2D
#define __BasketPut
//#define __BasketProdPut
//#define __SumOfPuts

using namespace std;
using namespace lawa;

typedef long double T;

typedef flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >  DenseMatrixT;
typedef flens::DenseVector<flens::Array<T> >                        DenseVectorT;

///  Basis definitions
typedef Basis<T,Orthogonal,Interval,Multi>                          PrimalBasis;
typedef TensorBasis2D<Adaptive,PrimalBasis,PrimalBasis>             Basis2D;

/* ************************************ */
/* *** Typedefs for financial model *** */
/* ************************************ */

///  Values for strike, maturity and weights for the assets $S_1$, $S_2$
T strike = 1.;
T maturity = 1;
T weight1 = 0.5, weight2 = 0.5;




///  Definition of the process type: Here (for the moment) the two-dimensional CGMYe model (see Section 8.7.1)


//T r = 0.04; T sigma1 = 0.3, sigma2 = 0.2, rho = 0.;
//T u11 = 1., u12 = 0., u21 = 0., u22 = 1.;

///  Parameters for the diffusion part
T r = 0.03; //0.03; //LK, was 0.03
T sigma1 = 0.3; //LK, was 0.3
T sigma2 = 0.3; //LK, was 0.3, was 0.2
T rho = 0.;

//T k_C1 = 1., k_G1 = 7.4, k_M1 = 8.5, k_Y1 = 0.8;
//T k_C2 = 1., k_G2 = 6.5, k_M2 = 9.5, k_Y2 = 1.1;

///  Parameters for the CGMY process parts
// LK: is this a martingale with the diffusion part? -> not necessarily, drift chosen s.t. martingale
//LK, wasT k_C1 = 1., k_G1 = 8.7, k_M1 = 16.5, k_Y1 = 1.25;
//LK, wasT k_C2 = 1., k_G2 = 11.2, k_M2 = 7.9, k_Y2 = 1.55;

//Lk
// for SumOfPutsSigma1111
//T k_C1 = 1., k_G1 = 16, k_M1 = 8, k_Y1 = 1.5; // 1
// -- end for basketprodput
//T k_C1 = 1., k_G1 = 8, k_M1 = 16, k_Y1 = 1.5; //1
T k_C1 = 1., k_G1 = 8, k_M1 = 16, k_Y1 = 0.2; //2
//T k_C1 = 1., k_G1 = 8, k_M1 = 8, k_Y1 = 1.5; //3
//T k_C1 = 1., k_G1 = 16, k_M1 = 16, k_Y1 = 1.5; //4
T k_C2 = 1., k_G2 = 16, k_M2 = 8, k_Y2 = 1.5; // 1
//T k_C2 = 1., k_G2 = 4, k_M2 = 4, k_Y2 = 1.5; // Analysis2

//T k_C2 = 1., k_G2 = 8, k_M2 = 8, k_Y2 = 1.25; // für 1-d vergleich
//T k_C2 = k_C1, k_G2 = k_G1, k_M2=k_M1, k_Y2 = k_Y1;

//LK: Transformation Matrix, max(||) <! G,M
//T Sigma11=.95, Sigma12=0.05, Sigma21=.05, Sigma22=.95;
//T Sigma11=1., Sigma12=1., Sigma21=1, Sigma22=1.;
//T Sigma11=.5, Sigma12=.4, Sigma21=.4, Sigma22=.5;
//T Sigma11=.5, Sigma12=.5, Sigma21=.5, Sigma22=.5;

//LK: Sigma_1
//T Sigma11=.9, Sigma12=0.1, Sigma21=0.2, Sigma22=.8;
//LK: Sigma_2
//T Sigma11=.7, Sigma12=0.3, Sigma21=0.4, Sigma22=.6;
//LK: Sigma_4
//T Sigma11=1., Sigma12=0., Sigma21=0., Sigma22=1.;

T alpha = 0.39;
T Sigma11 = 1-alpha, Sigma12 = alpha, Sigma21 = alpha, Sigma22 = 1-alpha;

//T Sigma11=1., Sigma12=0.2, Sigma21=0.3, Sigma22=1.;
//T Sigma11=.6, Sigma12=0.4, Sigma21=0.3, Sigma22=.7;
//T Sigma11=.95, Sigma12=0.05, Sigma21=0.05, Sigma22=.95;
//T Sigma11=.9, Sigma12=0.1, Sigma21=0.1, Sigma22=0.9;
//T Sigma11=.8, Sigma12=0.2, Sigma21=0.2, Sigma22=0.8;
//T k_C1 = 1., k_G1 = 8.7, k_M1 = 16.5, k_Y1 = 0.2;
//T k_C2 = 1., k_G2 = 11.2, k_M2 = 7.9, k_Y2 = 1.55;
//T k_C1 = 1., k_G1 = 8.7, k_M1 = 16.5, k_Y1 = 1.55;
//T k_C2 = 1., k_G2 = 11.2, k_M2 = 7.9, k_Y2 = 1.55;

// These parameters are required later on for computing wavelet integrals against the payoff function
// in order to determine when the payoff function is zero

//LK, since there is a linear transformation of the domain, there are no axes (only transformed axes) that can be used
// to determine whether the integral on a rectangle is zero (at least not with a minimum effort)
// improvement possible
//T    critical_line_x1 = 0.6;
//bool critical_above_x1 = false;//LK, was true;


//LK to store output folder
string outputFolder;


//LK  ----  process and option parameters ----------------------------------------------------------------------------
///  Storing the process parameters


#ifdef __CGMYeUnivariateJump2D
  ProcessParameters2D<T,CGMYeUnivariateJump2D>   processparameters(r, sigma1, sigma2, rho,
																	 k_C1,  k_G1, k_M1, k_Y1,
																	 k_C2,  k_G2, k_M2, k_Y2);
  const ProcessType2D  processtype  = CGMYeUnivariateJump2D;
#else
//#ifdef __CGMYeLinearTransformation2DJump2D
  ProcessParameters2D<T,CGMYeLinearTransformation2DJump2D>   processparameters(r, sigma1, sigma2, rho,
 																				 k_C1,  k_G1, k_M1, k_Y1,
 																				 k_C2,  k_G2, k_M2, k_Y2,
 																				 Sigma11, Sigma12,
 																				 Sigma21, Sigma22);
  const ProcessType2D  processtype  = CGMYeLinearTransformation2DJump2D;
#endif

///  Definition of the option type.
#ifdef __SumOfPuts
//#ifdef __CGMYeUnivariateJump2D
  // 	//LK only for independent underlyings!
 	OptionParameters2D<T,SumOfPuts> optionparameters(strike, strike, maturity, weight1, weight2, false);
// #endif
  const OptionTypenD optiontype = SumOfPuts;
  typedef PayoffIntegral2D<FullGridGL,Basis2D,TruncatedSumOfPutsOption2D<T> > PayoffIntegral;
#endif

#ifdef __BasketPut
  #ifdef __CGMYeLinearTransformation2DJump2D
    OptionParameters2D<T,BasketPut> optionparameters(strike, maturity, weight1, weight2,
 													 false, true, processparameters);   
  #endif
  #ifdef __CGMYeUnivariateJump2D
		OptionParameters2D<T,BasketPut> optionparameters(strike, maturity, weight1, weight2, false);
  #endif
  const OptionTypenD optiontype = BasketPut;
  typedef PayoffIntegral2D<FullGridGL,Basis2D,TruncatedBasketPutOption2D<T> > PayoffIntegral;
#endif

#ifdef __BasketProdPut
  #ifdef __CGMYeLinearTransformation2DJump2D
    OptionParameters2D<T,BasketProdPut> optionparameters(strike, maturity,false);   
  #endif
  #ifdef __CGMYeUnivariateJump2D
		OptionParameters2D<T,BasketProdPut> optionparameters(strike, maturity, false);
  #endif
  const OptionTypenD optiontype = BasketProdPut;
  typedef PayoffIntegral2D<FullGridGL,Basis2D,TruncatedBasketProdPutOption2D<T> > PayoffIntegral;
#endif
//LK ----  end process and option parameters ---------------------------------------------------------------------------

/* ********************************************* */
/* *** Typedefs for numerical discretization *** */
/* ********************************************* */

//typedef OptimizedH1Preconditioner2D<T,Basis2D>                      Preconditioner;

///  Definition of the underlying operator for the two-dimensional CGMYe process



#ifdef __CGMYeUnivariateJump2D
  typedef FinanceOperator2D<CGMYeUnivariateJump2D, Basis2D>           CGMYeOp2D;
#endif

#ifdef __CGMYeLinearTransformation2DJump2D
  typedef FinanceOperator2D<CGMYeLinearTransformation2DJump2D, Basis2D> CGMYeOp2D;
#endif

///  Local operator for the time-stepping scheme (see, e.g., Eq. (8.73))
typedef ThetaTimeStepLocalOperator<Index2D,CGMYeOp2D>               ThetaTimeStepLocalOperator2D;

///  Preconditioner adapted to the CGMYe operator
typedef DiagonalMatrixPreconditioner2D<T,Basis2D,CGMYeOp2D>         Preconditioner;

///  Required right-hand side definitions
typedef RHSWithPeaks1D<T,PrimalBasis>                               Rhs1D;
typedef AdaptiveSeparableRhs<T,Index2D,Rhs1D,Rhs1D >                AdaptiveSeparableRhsIntegral2D;
typedef ThetaTimeStepSeparableRHS<T,Index2D,
                                  AdaptiveSeparableRhsIntegral2D,
                                  ThetaTimeStepLocalOperator2D>     ThetaTimeStepRhs2d;
typedef CompoundRhs<T,Index2D,AdaptiveSeparableRhsIntegral2D,
                    AdaptiveSeparableRhsIntegral2D>                 CompoundRhsIntegral2D;

///  Definition of multitree AWGM solver for solving the linear system in each time-step. Here, we
///  we will only perform one AWGM step which will correspond to a sparse grid scheme. Note that
///  this proceeding can be optimized by writing a separate routine that does not rely on the
///  multitree solver class
typedef MultiTreeAWGM<Index2D,Basis2D,ThetaTimeStepLocalOperator2D,
                      ThetaTimeStepRhs2d,Preconditioner>            ThetaTimeStepMultiTreeAWGM2D;

///  Definition of multitree AWGM solver for approximating the initial condition. Here, we
///  we will only perform one AWGM step which will correspond to a sparse grid scheme. Note that
///  this proceeding can be optimized by writing a separate routine that does not rely on the
///  multitree solver class
typedef MultiTreeAWGM<Index2D,Basis2D,
                      ThetaTimeStepLocalOperator2D,
                      CompoundRhsIntegral2D,
                      NoPreconditioner<T,Index2D> >                 ApproxL2AWGM2D;

///  Definition of the $\theta$-scheme AWGM sovler
typedef ThetaSchemeAWGM<Index2D, ThetaTimeStepMultiTreeAWGM2D>      ThetaSchemeMultiTreeAWGM2D;



///  Iterators for post-processing
typedef IndexSet<Index1D>::const_iterator                           const_set1d_it;
typedef IndexSet<Index2D>::const_iterator                           const_set2d_it;
typedef Coefficients<Lexicographical,T,Index1D>::iterator           coeff1d_it;
typedef Coefficients<Lexicographical,T,Index2D>::iterator           coeff2d_it;
typedef Coefficients<Lexicographical,T,Index1D>::const_iterator     const_coeff1d_it;
typedef Coefficients<Lexicographical,T,Index2D>::const_iterator     const_coeff2d_it;


T f_t(T t)       {  return 0.; }
T f_x(T x)       {  return 0.; }
T f_y(T y)       {  return 0.; }


/// A simple routine to evaluate a wavelet basis expansion on the domain $[-R_1,R_1] \times [-R_2,R_2]$
T
evaluate(const Basis2D &basis2d, T left_x1, T right_x1, T left_x2, T right_x2,
         const Coefficients<Lexicographical,T,Index2D> &v, T x1, T x2);

///  Computing the $L_\infty$ error as described in Eq. (8.125). There, you also find the definition
///  of $\delta$.



T
computeLinftyError(const Basis2D &basis2d, T left_x1, T right_x1, T left_x2, T right_x2,
                   const Coefficients<Lexicographical,T,Index2D> &u,T delta, int j,
                   Option2D<T,optiontype> &option2d,
                   ProcessParameters2D<T,processtype> &processparameters);



///  For scatter plots (=patterns) of the stiffness matrix
void
spyStiffnessMatrix(const Basis2D &basis2d, CGMYeOp2D & cgmyeop2d, int j,
                   const ProcessParameters2D<T,processtype> &processparameters);


void
computeSensitivitiesDelta1and2(const Basis2D &basis2d, T left_x1, T right_x1, T left_x2, T right_x2,
							   const Coefficients<Lexicographical,T,Index2D> &u, const Coefficients<Lexicographical,T,Index2D> &u_nextToLast, T dt,
							   T delta, int j, Option2D<T,optiontype> &option2d,
							   ProcessParameters2D<T,processtype> &processparameters);



int main (int argc, char *argv[]) {

    cout.precision(16);
    if ((argc<=11) || (argc > 13)) {
        cout << "Usage: " << argv[0] << " d j0 J R1_1 R2_1 R1_2 R2_2 N(MC) m(time) delta(error) compSensitivities folder(Output,optional)" << endl;
        return 0;
    }

    ///  Wavelet basis parameters
    int d                = atoi(argv[1]);
    int j0               = atoi(argv[2]);
    int J                = atoi(argv[3]);
	int N                = atoi(argv[8]);
	int numOfTimesteps   = atoi(argv[9]);
	double delta         = atof(argv[10]);
	int compSensint      = atoi(argv[11]);
	outputFolder= "";
	if (argc==13){
		outputFolder=argv[12];
		if (outputFolder[outputFolder.length()-1] != '/') outputFolder= outputFolder +'/';//strcat(outputFolder,"/");
	}
	
	bool computeSensitivities = (compSensint == 1) ? true : false; 

    ///  Parameters for AWGM: here only dummy variables
    T alpha = 0.7;
    T gamma = 0.025;
    const char* residualType = "standard";
    const char* treeType = "sparsetree"; //"gradedtree";

    ///  We focus on $L_2$-orthonormal wavelets here
    bool IsMW = true;


    // LK: why fix?? Guess: don't know size before, so take very large one
    size_t hashMapSize = 196613;

    ///  Size of the underlying domain: $[-R_1,R_1] \times [-R_2,R_2]$
    T R1_1 = atof(argv[4]);
    T R2_1 = atof(argv[5]);
    T left_x1 = -R1_1, right_x1 = R2_1;
    T R1_2 = atof(argv[6]);
    T R2_2 = atof(argv[7]);
    T left_x2 = -R1_2, right_x2 = R2_2;

    ///  Parameter $\delta$ for error measurement (e.g., Eq. (8.125))
    //LK: was T delta = 0.05;
    //T delta = 1.0;

    ///  Parameters for the $\theta$-scheme
    T theta = 0.5;
    T timestep_eps = 1e-6;
    int maxiterations =  1;  T init_cgtol = 1e-9;   // use maxiterations = 1 for "pure" sparse grid computation
    //int numOfTimesteps = 128;
    T timestep = maturity/numOfTimesteps;

    ///  Integration order for approximating the initial condition
    int order = 4;

    ///  Read reference prices from file (true) or not (false)
    bool useRefPrices = false;

    Timer time;
	
    /// Basis initialization
    PrimalBasis     basis(d,j0);
    basis.enforceBoundaryCondition<DirichletBC>();
    Basis2D         basis2d(basis,basis);

	///  Initialization of the CGMY operator. Observe that we need to pass the parameters of the
	///  underlying domain as parameters for a domain transformation!
	CGMYeOp2D                    cgmyeOp2D(basis2d, processparameters,
										   R1_1, R2_1, R1_2, R2_2, 10);
	std::cout << "Processparameters = " << processparameters << std::endl;
	///  Initialization of the time-step oerator
	ThetaTimeStepLocalOperator2D localThetaTimeStepOp2D(theta,timestep,cgmyeOp2D);
    
	/// Initialization of preconditioner
	NoPreconditioner<T, Index2D> NoPrec;
	//Preconditioner  Prec(basis2d, sigma1*sigma1, sigma2*sigma2, 1.);
	Preconditioner  Prec(cgmyeOp2D);


	// Initialization of integrals for the rhs which is, for our example, zero. The example below
	// however shows how to use singular points for a refinement of the integration domain when
	// the function to be integrated against is not smooth at or near the origin. Please note that
	// such a rhs object is required for the implementation of the $\theta$-scheme which also applies
	// to more general problems
	DenseVectorT sing_pts_t, sing_pts_x, sing_pts_y;
	//sing_pts_x = 0.1, 0.2, 0.3, 0.4, 0.5;
	//sing_pts_y =  0.1, 0.2, 0.3, 0.4, 0.5;
	DenseMatrixT no_deltas, deltas_x, deltas_y;
	Function<T>                    fct_f_t(f_t,sing_pts_t);
	Function<T>                    fct_f_x(f_x,sing_pts_x), fct_f_y(f_y,sing_pts_y);
	RHSWithPeaks1D<T,PrimalBasis>  rhs_f_x(basis, fct_f_x, no_deltas, order);
	RHSWithPeaks1D<T,PrimalBasis>  rhs_f_y(basis, fct_f_y, no_deltas, order);
	Coefficients<Lexicographical,T,Index1D> rhs_f_x_data(SIZEHASHINDEX1D),
		rhs_f_y_data(SIZEHASHINDEX1D);
	AdaptiveSeparableRhsIntegral2D F_rhs(rhs_f_x, rhs_f_x_data, rhs_f_y, rhs_f_y_data);
	ThetaTimeStepRhs2d thetatimestep_F(fct_f_t,F_rhs,localThetaTimeStepOp2D);
	// LK: how to do thetatimestep_F(...)??


    
	/// Initialization of integrals for initial condition angd rhs
	Option2D<T,optiontype>         option2d(optionparameters);

	//LK: set linear transformation
    #ifdef  __CGMYeLinearTransformation2DJump2D//__BasketProdPut
	  #ifdef __BasketPut
	    option2d.calc_w();
	  #else
	    option2d.calc_w(processparameters);
	  #endif
    #endif
	
	//LK: set MC-NrOfSimulations
	#ifdef __BasketPut
        option2d.setLevyMeasureTruncationEps(0.000001); //(0.000001)
		option2d.setNumberOfMCRuns(N); // was 100000
    #endif


	///  This is required for approximating the initial condition with zero boundary conditions (see, e.g., p. 183)
    #ifdef __BasketPut
	    TruncatedBasketPutOption2D<T> truncatedoption2d;
    #endif

    #ifdef __SumOfPuts
	    TruncatedSumOfPutsOption2D<T> truncatedoption2d;
    #endif

	#ifdef __BasketProdPut
	    TruncatedBasketProdPutOption2D<T> truncatedoption2d;
    #endif	
    	

	// LK: check ---------------------------------------------------------------------------------------------------------------
     
	truncatedoption2d.setOption(option2d);
	// LK: (left_xy,right_xy, left_x2,right_x2, _type, _truncWidth, damping_c)
	truncatedoption2d.setTruncation(left_x1, right_x1, left_x2, right_x2, 0, 0.1, 100.);
	//truncatedoption2d.setCriticalLine_x1(critical_line_x1, critical_above_x1);
	PayoffIntegral payoffIntegral(basis2d, truncatedoption2d,
								  left_x1, right_x1, left_x2, right_x2, true, 0.05, order);
	//LK_------------------------------------------------------------------------------------------------------------------------


	
	//LK: initializing u,f
	Coefficients<Lexicographical,T,Index2D> u(hashMapSize), f(hashMapSize), u_nextToLast(hashMapSize);
    
	std::stringstream filename;

	if (optiontype == BasketPut) {
		filename << outputFolder << "LinTrans_bput2d_conv2_sg_" << d << "_" << "_lx1_" << left_x1 << "_rx1_" << right_x1
				 << "_lx2_" << left_x2 << "_rx2_" << right_x2 << "_del_" << delta 
				 << "_strike_" << strike << "_maturity_" << maturity << "_weight1_" << weight1 << "_weight2_" << weight2
				 << processparameters << ".txt";
	}
	else if (optiontype == SumOfPuts) {
		filename << outputFolder << "LinTrans_sumofput2d_sg_" << d << "_" << "_lx1_" << left_x1 << "_rx1_" << right_x1
				 << "_lx2_" << left_x2 << "_rx2_" << right_x2 << "_delta_" << delta
				 << "_strike_" << strike << "_maturity_" << maturity << "_weight1_" << weight1 << "_weight2_" << weight2
				 << processparameters << ".txt";
	}

	std::ofstream convfile(filename.str().c_str());
		
	for (int j=0; j<=J; ++j) {
		time.start();
		
		getSparseGridVector(basis2d, u, j, (T)0.);

		///  This required for the correct initialization of the compression for CGMY operator (see p. 186)
		cgmyeOp2D.setCompressionLevel(j, j);

		//spyStiffnessMatrix(basis2d, cgmyeOp2D, j, processparameters);
		cerr << "Computation of initial condition started." << endl;
		
		int count = 0;

		for (coeff2d_it it=u.begin(); it!=u.end(); ++it) {
			coeff2d_it it_f = f.find((*it).first); //LK: f used for recalculating (next j)
			if (it_f != f.end()) { //LK: first time -> always else-case
				(*it).second = (*it_f).second;
			}
			else {
			    
				T tmp = payoffIntegral((*it).first);
				f[(*it).first] = tmp;
				(*it).second = tmp;
			}
			++count;
			if ( (count%(10* ((int)std::pow(2,j) ))==1) || (count==u.size()) ) cout << "count: " << count << " / " << u.size() << endl;
		}
		timestep = maturity/numOfTimesteps;

		/// Initialization of multi tree based adaptive wavelet Galerkin method
		ThetaTimeStepMultiTreeAWGM2D thetatimestep_solver(basis2d, localThetaTimeStepOp2D,
														  thetatimestep_F, Prec);
		thetatimestep_solver.setParameters(alpha, gamma, residualType, treeType, IsMW, false,
										   hashMapSize);

		/// Initialization of $\theta$-scheme solver. The value of "zero" (last argument) refers
		/// to a sparse grid like realization of AWGM used in each time-step. More precisely, if
		/// we set this option, the initial index set associated to the wavelet basis expansion of
		/// the initial condition is fixed and not changed in the course of the time-stepping.
		ThetaSchemeMultiTreeAWGM2D thetascheme(thetatimestep_solver); 
		thetascheme.setParameters(theta, timestep, numOfTimesteps, timestep_eps, maxiterations,
								  init_cgtol, 0); // LK: strategy = 0 // LK, true= saves forelast solution (for sensitivites wrt t)

		///  Dummy variables: these are only needed for adaptive computations. In this program,
		///  we only use sparse grid index sets
		int avDof = 0, maxDof = 0., terminalDof;

		///  Calling the $\theta$-scheme solver
		//LK, was thetascheme.solve(u, avDof, maxDof, terminalDof, j); 
		thetascheme.solveWithNextToLast(u, u_nextToLast, avDof, maxDof, terminalDof, j);  // LK -> adapt so that u_forelast saves the forlast solution saves u(T-dt)
				
		cerr << "Computation of u has finished." << endl;
		T maxerror = 0., maxerror1 = 0., maxerror2 = 0.;
		maxerror = computeLinftyError(basis2d, left_x1, right_x1, left_x2, right_x2, u, delta, j, option2d,
									  processparameters);

		if (computeSensitivities){
		computeSensitivitiesDelta1and2(basis2d, left_x1, right_x1, left_x2, right_x2, u, u_nextToLast,
									   timestep, delta, j, option2d, processparameters);
		}
		//spyStiffnessMatrix(basis2d, cgmyeOp2D, j, processparameters);
		
		convfile << timestep << " " << j << " " << u.size() << " "
			//                                          theoretical conv       MC conv order
     		     << maxerror << " " << delta <<" " << std::pow(4,-j) << " " << 1./std::sqrt(N)  <<endl;
		cerr << "Computation of errors has finished." << endl;
		//computeReferencePrice(basis2d, left_x1, right_x1, left_x2, right_x2,
		//                      -0.1, 0.1, -0.1, 0.1, 0.02, 0.02, u, j, option2d, processparameters);
		cerr << "Computation of reference prices has finished." << endl;
		// LK, solution V(1,1)
		//std::cout << " Main-programm, w1 = " << option2d.get_w1() << ", w2 = " << option2d.get_w2() << std::endl;
		T x1hat = 0.0;
		T x2hat =0.0;
		#ifdef __CGMYeUnivariateJump2D
		std::cout << "univariate jumps" << std::endl;
		    x1hat = r*maturity;
		    x2hat = r*maturity;
        #endif
			T approx = std::exp(-r*maturity)* evaluate(basis2d, left_x1, right_x1, left_x2, right_x2, u, x1hat, x2hat);
		cout << "j = " << setw(2) << j << setw(18) << left << ", Solution V(1,1,0) = " << setw(18)<<  approx ;
		T exact = 0.0;

        #ifdef __BasketProdPut 
		  exact = option2d.value(processparameters,1.,1.,0);
		#else
		  exact = option2d.value(processparameters,1.,1.,0);
		#endif
		time.stop();
		cout << ", MC-Solution V(1,1,0) = " << setw(18) << left <<  exact
			 <<", duration of total calculation: " << std::floor(time.elapsed()/((double)60))
			 << " m, " << ((int)std::floor(time.elapsed())) % 60 << " sec"<< std::endl;
				
	}

    return 0;
}



T
evaluate(const Basis2D &basis2d, T left_x1, T right_x1, T left_x2, T right_x2,
         const Coefficients<Lexicographical,T,Index2D> &v, T x1, T x2)
{
    T RightmLeft_x1 = right_x1-left_x1, SqrtRightmLeft_x1 = std::sqrt(right_x1-left_x1);
    T RightmLeft_x2 = right_x2-left_x2, SqrtRightmLeft_x2 = std::sqrt(right_x2-left_x2);

    T ret = 0.;

    for (const_coeff2d_it it=v.begin(); it!=v.end(); ++it) {
        int   j1 = (*it).first.index1.j,     j2 = (*it).first.index2.j;
        int   k1 = (*it).first.index1.k,     k2 = (*it).first.index2.k;
        XType e1 = (*it).first.index1.xtype, e2 = (*it).first.index2.xtype;

        T val_x1 = (1./SqrtRightmLeft_x1) * basis2d.first.generator(e1).operator()((x1-left_x1)/(RightmLeft_x1),j1,k1,0);
        T val_x2 = (1./SqrtRightmLeft_x2) * basis2d.second.generator(e2).operator()((x2-left_x2)/(RightmLeft_x2),j2,k2,0);

        ret += (*it).second * val_x1 * val_x2;
    }
    return ret;
}






T
computeLinftyError(const Basis2D &basis2d, T left_x1, T right_x1, T left_x2, T right_x2,
                   const Coefficients<Lexicographical,T,Index2D> &u,T delta, int j,
                   Option2D<T,optiontype> &option2d,
                   ProcessParameters2D<T,processtype> &processparameters)
{
	Timer time;
	
    std::stringstream filename;
    if (optiontype == BasketPut) {
        filename //<< "cgmye2d_basketput_sg_" << j << "_lx1_" << left_x1 << "_rx1_" << right_x1
			<< outputFolder << "cgmye2d_bput_sg_" << j << "_lx1_" << left_x1 << "_rx1_" << right_x1
                 << "_lx2_" << left_x2 << "_rx2_" << right_x2 << "_delta_" << delta
                 << "_strike_" << strike << "_maturity_" << maturity << "_weight1_" << weight1 << "_weight2_" << weight2
                 << processparameters << ".txt";
    }
    else if (optiontype == SumOfPuts) {
        filename << outputFolder << "cgmye2d_sumputs_sg_" << j << "_lx1_" << left_x1 << "_rx1_" << right_x1
                 << "_lx2_" << left_x2 << "_rx2_" << right_x2 << "_delta_" << delta
                 << "_strike_" << strike << "_maturity_" << maturity << "_weight1_" << weight1 << "_weight2_" << weight2
                 << processparameters << ".txt";
    }
    else {
        std::cerr << "Unknown option type" << std::endl;
		exit(1);
    }
    std::ofstream plotfile(filename.str().c_str());
    plotfile.precision(16);

	T maxerror = 0.;
    /*
	  T _left_x1  = left_x1, _left_x2  = left_x2;
	  T _right_x1 = left_x1, _right_x2 = left_x2;
    */

    T _left_x1  = std::min(left_x1,left_x2);
    T _left_x2  = _left_x1;
    T _right_x1 = std::max(right_x1,right_x2);
    T _right_x2 = _right_x1;

    T h1 = (delta*_right_x1-delta*_left_x1)/32.;
	T h2 = (delta*_right_x1-delta*_left_x1)/32;
	

    std::cerr << "ATTENTION: adapted error measuremnt!!" << std::endl;

    std::cerr << "Error is measured over [" << delta*_left_x1 << ", " << delta*_right_x1 << "], ["
              << delta*_left_x2 << ", " << delta*_right_x2 << "]" << std::endl;
    //for (T x1=left_x1; x1<=right_x1; x1+=0.03125) {
    for (T x1=delta*_left_x1; x1<=delta*_right_x1; x1+=h1) {
		
        for (T x2=delta*_left_x2; x2<=delta*_right_x2; x2+=h2) {
			
            T S1    = strike*std::exp(x1);
            T S2    = strike*std::exp(x2);

			T Sigmax1 = Sigma11 * x1 + Sigma12 * x2;
			T Sigmax2 = Sigma21 * x1 + Sigma22* x2;

			T SSigmax1 = std::exp(Sigmax1);
			T SSigmax2 = std::exp(Sigmax2);
			time.start();
			//S1,S2 = exp(Y_0)
			T exact=0.0;
			#ifdef __BasketProdPut
			exact = option2d.value(processparameters,S1,S2,0);
			#else
			  exact = option2d.value(processparameters,S1,S2,0);
			#endif
			time.stop();
			if ((x1==delta*_left_x1) && (x2==delta* _left_x2)) std::cout <<" reference price calculation duation is " << time.elapsed() << " seconds " << std::endl;
			T x1hat = x1;
			T x2hat = x2;
            						
			//payoff (now (S_0))
            T payoff = 0.0;

			#ifdef __CGMYeLinearTransformation2DJump2D
			  #ifdef __BasketPut
			    payoff = option2d.payoff(strike*exp(x1),strike*exp(x2), - (r + option2d.get_w1())*maturity,-(r+option2d.get_w2())*maturity);
				//#endif
			  #else
			    payoff = option2d.payoff(strike*exp(x1),strike*exp(x2),(T) 0.);
				//payoff = option2d.payoff_log(x1,x2);
				//std::cout << " exp(x1) = " << std::exp(x1) << " exp(x2) = " << std::exp(x2) << "payoff = " << payoff  << std::endl;
			  #endif			  
            #endif
			  //LK schön machen!

			#ifdef __CGMYeUnivariateJump2D  
			  payoff = option2d.payoff(strike*exp(x1),strike*exp(x2)); //LK check
			#endif

			  
			if (processtype==CGMYeUnivariateJump2D){
				x1hat+= r*maturity; //LK was - 
				x2hat+= r*maturity; //LK was -
			}

			//LK: anpassen!
			// x1hat,x2hat beziehen sich auf N+Y 
            T approx =std::exp(-r*maturity)*evaluate(basis2d, left_x1, right_x1, left_x2, right_x2, u, x1hat, x2hat);

			//                         1                         2                        3               4               5              6                7
			plotfile	<< std::exp(Sigmax1) << " " << std::exp(Sigmax2) << " " << Sigmax1 << " " << Sigmax2 << " " << exact << " " << approx << " " << payoff 
				//      8                            9  
				//" "  << Sigmax1hat << " " << Sigmax2hat	<<
						<< " " << std::endl;
            maxerror = std::max(maxerror, fabs(approx-exact));
        }
        plotfile << endl;
    }
    plotfile.close();

    return maxerror;
}

void
spyStiffnessMatrix(const Basis2D &basis2d, CGMYeOp2D & cgmyeop2d, int j,
                   const ProcessParameters2D<T,processtype> &processparameters)
{
    typedef flens::SparseGeMatrix<flens::CRS<T,flens::CRS_General> >  SparseMatrixT;
    typedef std::list<Index2D>::const_iterator                        const_list_it;
    std::list<Index2D> row_indices, col_indices, indices;

    int j0_x = basis2d.first.j0;
    int j0_y = basis2d.second.j0;

    for (int levelsum=0; levelsum<=j; ++levelsum) {
        for (int i1=0; i1<=levelsum; ++i1) {
            for (int i2=0; i2<=levelsum; ++i2) {
                if (i1+i2!=levelsum) continue;

                if ((i1==0) && (i2==0)) {
                    //cout << "(" << i1 << ", " << i2 << ") : " << basis2d.first.mra.cardI(j0_x) * basis2d.second.mra.cardI(j0_y) << endl;
                    for (long k1=basis2d.first.mra.rangeI(j0_x).firstIndex(); k1<=basis2d.first.mra.rangeI(j0_x).lastIndex(); ++k1) {
                       for (long k2=basis2d.second.mra.rangeI(j0_y).firstIndex(); k2<=basis2d.second.mra.rangeI(j0_y).lastIndex(); ++k2) {
                           Index1D row(j0_x,k1,XBSpline);
                           Index1D col(j0_y,k2,XBSpline);
                           indices.push_back(Index2D(row,col));
                       }
                    }
                }
                else if ((i1!=0) && (i2==0)) {
                    int j1=j0_x+i1-1;
                    //cout << "(" << i1 << ", " << i2 << ") : " << basis2d.first.cardJ(j1) * basis2d.second.mra.cardI(j0_y) << endl;
                    for (long k1=basis2d.first.rangeJ(j1).firstIndex(); k1<=basis2d.first.rangeJ(j1).lastIndex(); ++k1) {
                        for (long k2=basis2d.second.mra.rangeI(j0_y).firstIndex(); k2<=basis2d.second.mra.rangeI(j0_y).lastIndex(); ++k2) {
                            Index1D row(j1,k1,XWavelet);
                            Index1D col(j0_y,k2,XBSpline);
                            indices.push_back(Index2D(row,col));
                        }
                    }
                }
                else if ((i1==0) && (i2!=0)) {
                    int j2=j0_y+i2-1;
                    //cout << "(" << i1 << ", " << i2 << ") : " << basis2d.first.mra.cardI(j0_x) * basis2d.second.cardJ(j2) << endl;
                    for (long k1=basis2d.first.mra.rangeI(j0_x).firstIndex(); k1<=basis2d.first.mra.rangeI(j0_x).lastIndex(); ++k1) {
                        for (long k2=basis2d.second.rangeJ(j2).firstIndex(); k2<=basis2d.second.rangeJ(j2).lastIndex(); ++k2) {
                            Index1D row(j0_x,k1,XBSpline);
                            Index1D col(j2,k2,XWavelet);
                            indices.push_back(Index2D(row,col));
                        }
                    }
                }
                else if ((i1!=0) && (i2!=0)) {
                    int j1=j0_x+i1-1;
                    int j2=j0_y+i2-1;
                    //cout << "(" << i1 << ", " << i2 << ") : " << basis2d.first.cardJ(j1) * basis2d.second.cardJ(j2) << endl;
                    for (long k1=basis2d.first.rangeJ(j1).firstIndex(); k1<=basis2d.first.rangeJ(j1).lastIndex(); ++k1) {
                        for (long k2=basis2d.second.rangeJ(j2).firstIndex(); k2<=basis2d.second.rangeJ(j2).lastIndex(); ++k2) {
                            Index1D row(j1,k1,XWavelet);
                            Index1D col(j2,k2,XWavelet);
                            indices.push_back(Index2D(row,col));
                        }
                    }
                }
            }
        }
    }


    int N = indices.size();
    //std::cerr << "   N = " << N << std::endl;
    SparseMatrixT A(N,N);

    int row_pos = 1, col_pos = 1;
    for (const_list_it row_it=indices.begin(); row_it!=indices.end(); ++row_it) {
        for (const_list_it col_it=indices.begin(); col_it!=indices.end(); ++col_it) {
            T val = cgmyeop2d((*row_it),(*col_it));
            if (fabs(val)>0) A(row_pos,col_pos) = val;
            ++col_pos;
        }
        ++row_pos;
        col_pos = 1;
    }
    cout << "A finalize called..." << endl;
    A.finalize();
    cout << "... finished" << endl;
    std::stringstream matrixfilename;
    matrixfilename << outputFolder  << "A_" << N << "_" << processparameters;
    cout << "Spy called..." << endl;
    spy(A,matrixfilename.str().c_str(),true,(T)0.);
    cout << "... finished" << endl;
}





//LK, computing Delta^i= d/dS_i V(t,S) and Theta = d/dt V(t,S)
void
computeSensitivitiesDelta1and2(const Basis2D &basis2d, T left_x1, T right_x1, T left_x2, T right_x2,
							   const Coefficients<Lexicographical,T,Index2D> &u,
							   const Coefficients<Lexicographical,T,Index2D> &u_nextToLast, T dt,
							   T delta, int j,
							   Option2D<T,optiontype> &option2d,
							   ProcessParameters2D<T,processtype> &processparameters)
{
	//computes delta S1,S2 at time zero
	//for different values S1,S2, such that $y1,y2 \in [left_x1*delta,right_x1*delta] \times [left_x2 * delta, right_x2 * delta]$
	Timer time;
		
    std::stringstream filename1;
	std::stringstream filename2;
	std::stringstream filename3;
	std::stringstream filename1_ref;
    if (optiontype == BasketPut) {
        filename1 //<< "cgmye2d_basketput_sg_" << j << "_lx1_" << left_x1 << "_rx1_" << right_x1
			<< outputFolder << "sensitivities_delta" << delta << "_DeltaS1_" << j << ".txt";
		//filename1_ref <<	outputFolder << "ref_sensitivities_delta_x1_" << j << ".txt";
		filename2 << outputFolder << "sensitivities_delta" << delta << "_DeltaS2_" << j << ".txt";
		filename3 << outputFolder << "sensitivities_delta" << delta << "_Theta_" << j << ".txt";
	}
    else {
        std::cerr << "Unknown option type for calculating sensitivities" << std::endl;
		return; //exit(1);
    }
    std::ofstream plotfile1(filename1.str().c_str());
	std::ofstream plotfile2(filename2.str().c_str());
	std::ofstream plotfile3(filename3.str().c_str());
	//std::ofstream plotfile1_ref(filename1_ref.str().c_str());

	
    plotfile1.precision(16);
	plotfile2.precision(16);
	plotfile3.precision(16);

    T _left_x1  = std::min(left_x1,left_x2);
    T _left_x2  = _left_x1;
    T _right_x1 = std::max(right_x1,right_x2);
    T _right_x2 = _right_x1;
	
	
    T h1 = (delta*_right_x1-delta*_left_x1)/32.;//((double)(std::pow(2,j+2)));
	T h2 = (delta*_right_x2-delta*_left_x2)/32.;//((double)(std::pow(2,j+2)));  // /32.;
	//T h = 0.0001;
	T h1_dq = (_right_x1-_left_x1)/std::pow(2,j); //h1 for differential quotient
	T h2_dq = (_right_x1-_left_x1)/std::pow(2,j); //h2 for differential quotient

	long step=1;
	long step_size = ((delta*_right_x1-delta*_left_x1)/h1) * ((delta*_right_x2-delta*_left_x2)/h2);

    std::cerr << "Sensitivities are approximated over [" << delta*_left_x1 << ", " << delta*_right_x1 << "], ["
              << delta*_left_x2 << ", " << delta*_right_x2 << "]" << std::endl;
	
	time.start();
	
    for (T x1=delta*_left_x1; x1<=delta*_right_x1; x1+=h1) {
        for (T x2=delta*_left_x2; x2<=delta*_right_x2; x2+=h2) {
			//if ((step++ % 1000) == 0)std::cout << " delta1 sensitivities step " << step  << " / " << step_size<< std::endl;
			T x1hat = x1;
			T x2hat = x2;
            
			if (processtype==CGMYeUnivariateJump2D){
				x1-= r*maturity;
				x2-= r*maturity;
			}

			//naiver weg
			T t=0;
			// u(x+h) - u(x-h)
			// -------------
			//      2h
			T u1 = evaluate(basis2d, left_x1, right_x1, left_x2, right_x2, u, x1hat - h1_dq, x2hat);
			T u2 = evaluate(basis2d, left_x1, right_x1, left_x2, right_x2, u, x1hat + h1_dq, x2hat);
			T u_y1 = (u2 - u1)/(2.*h1_dq);

			u1 = evaluate(basis2d, left_x1, right_x1, left_x2, right_x2, u, x1hat , x2hat - h2_dq);
			u2 = evaluate(basis2d, left_x1, right_x1, left_x2, right_x2, u, x1hat , x2hat + h2_dq);
			T u_y2 = (u2 - u1)/(2.*h2_dq);
			
            // Sigmax
			T Sigmax1 = Sigma11 * x1 + Sigma12 * x2;
			T Sigmax2 = Sigma21 * x1 + Sigma22 * x2;

			T S1 = std::exp(Sigmax1);
			T S2 = std::exp(Sigmax2);

			
			// Sigma has to be invertible!
			T det =  (Sigma11*Sigma22 - Sigma21*Sigma12);
			T V_S1 = std::exp(- r* (maturity -t)) * 1./(S1 * det)  * (Sigma22 * u_y1 - Sigma21 * u_y2);
			T V_S2 = std::exp(- r* (maturity -t)) * 1./(S2 * det)  *  (Sigma11 * u_y2 - Sigma12 * u_y1);
			/*
			// ---------------------- test x1 -----------------------------------
			// Sigma^-1 log(S) = y
			//T h1_ref = 0.0001;
			T h_ref= 0.00001;
			
			T S1ph = S1+h_ref;
			T S1mh = S1-h_ref;
			//T y1,y2;
			//T det = (Sigma11 * Sigma22 - Sigma21*Sigma12);

			// Sigma^-1 log(S)
			T y1_ref = 1./det * (Sigma22* std::log(S1) - Sigma12 * std::log(S2));
			T y2_ref = 1./det * (Sigma11* std::log(S2) - Sigma21 * std::log(S1));
			
			// Sigma^-1 (S1ph S2)^T
			T y1ph = 1./det * (Sigma22* std::log(S1ph) - Sigma12 * Sigmax2);
			T y2ph = 1./det * (Sigma11* Sigmax2 - Sigma21 * std::log(S1ph));

			T y1mh = 1./det * (Sigma22* std::log(S1mh) - Sigma12 * Sigmax2);
			T y2mh = 1./det * (Sigma11* Sigmax2 - Sigma21 * std::log(S1mh));


			T u1_ref = std::exp(-r *(maturity -t))* evaluate(basis2d, left_x1, right_x1, left_x2, right_x2, u, y1mh , y2mh);
			T u_ref =std::exp(-r *(maturity -t))* evaluate(basis2d, left_x1, right_x1, left_x2, right_x2, u, x1 , x2);
			T u2_ref =std::exp(-r *(maturity -t))*  evaluate(basis2d, left_x1, right_x1, left_x2, right_x2, u, y1ph , y2ph);
						
			T V_S1_ref = (u2_ref - u1_ref)/(2*h_ref); */

						
			// ---------------------- test x2-----------------------------------
			// Sigma^-1 log(S) = y
			//T h1_ref = 0.0001;
			/*
			T h_ref= 0.00001;
			
			T S2ph = S2+h_ref;
			T S2mh = S2-h_ref;
			//T y1,y2;
			//T det = (Sigma11 * Sigma22 - Sigma21*Sigma12);

			// Sigma^-1 log(S)
			T y1_ref = 1./det * (Sigma22* std::log(S1) - Sigma12 * std::log(S2));
			T y2_ref = 1./det * (Sigma11* std::log(S2) - Sigma21 * std::log(S1));
			
			// Sigma^-1 (S1 S2ph)^T
			T y1S2ph = 1./det * (Sigma22* std::log(S1) - Sigma12 * std::log(S2ph));
			T y2S2ph = 1./det * (Sigma11* std::log(S2ph) - Sigma21 * std::log(S1));
			// Sigma^-1 (S1 S2mh)^T
			T y1S2mh = 1./det * (Sigma22* std::log(S1) - Sigma12 * std::log(S2mh));
			T y2S2mh = 1./det * (Sigma11* std::log(S2mh) - Sigma21 * std::log(S1));


			
			T u1S2_ref = std::exp(-r *(maturity -t))* evaluate(basis2d, left_x1, right_x1, left_x2, right_x2, u, y1S2mh , y2S2mh);
			T u_ref =std::exp(-r *(maturity -t))* evaluate(basis2d, left_x1, right_x1, left_x2, right_x2, u, x1 , x2);
			T u2S2_ref =std::exp(-r *(maturity -t))*  evaluate(basis2d, left_x1, right_x1, left_x2, right_x2, u, y1S2ph , y2S2ph);
						
			T V_S2_ref = (u2S2_ref - u1S2_ref)/(2*h_ref); */
			
			// ------------- Theta -----------------------------------------------------
			// ----------- Method of choice ----------------;

			T u_Ty = evaluate(basis2d, left_x1, right_x1, left_x2, right_x2, u, x1 , x2);
			T u_Tmdty =  evaluate(basis2d, left_x1, right_x1, left_x2, right_x2, u_nextToLast, x1 , x2);

			T u_t = (-u_Tmdty + u_Ty)/dt;

			T V_t = 0.0;
			#ifdef __CGMYeLinearTransformation2DJump2D
			  #ifdef __BasketPut
				V_t = std::exp(-r * maturity) * ( r* u_Ty - u_t )  - (r + option2d.get_w1()) * S1 * V_S1 - (r+option2d.get_w2()) * S2 * V_S2;
			  #endif
			  #ifdef __BasketProdPut
			    V_t = std::exp(-r * maturity) * ( r* u_Ty - u_t )  - (option2d.get_wr1()) * S1 * V_S1 - (option2d.get_wr2()) * S2 * V_S2;
			  #endif			  
            #endif
			

			// ------------------ validation --------------
			// log(S) = Sigmax?
				/*
				T V_0S = std::exp(- r * maturity) * evaluate(basis2d, left_x1, right_x1, left_x2, right_x2, u, x1 , x2);

			// Sigma^-1 log(S)-(r+w)dt
			T logSmrpwtdt1 = 0.0;
			T logSmrpwtdt2 = 0.0;
			#ifdef __CGMYeLinearTransformation2DJump2D
			  #ifdef __BasketPut
			    logSmrpwtdt1 = Sigmax1 - (r+option2d.get_w1())*dt; //LK, was 0.0??
			    logSmrpwtdt2 = Sigmax2 - (r+option2d.get_w2())*dt;
			  #endif
			  #ifdef __BasketProdPut
				logSmrpwtdt1 = Sigmax1 - (option2d.get_wr1())*dt; //LK, was 0.0??
			    logSmrpwtdt2 = Sigmax2 - (option2d.get_wr2())*dt;
			  #endif			  
            #endif
			#ifdef __CGMYeLinearTransformation2DJump2D
				// LK, todo?
			#endif	

			T Sigmam1dt1 = 1./det * (Sigma22* logSmrpwtdt1 - Sigma12 * logSmrpwtdt2);
			T Sigmam1dt2 = 1./det * (Sigma11* logSmrpwtdt2 - Sigma21 * logSmrpwtdt1);

			T V_dtS = std::exp(- r* (maturity-dt)) * evaluate(basis2d, left_x1, right_x1, left_x2, right_x2, u_nextToLast, Sigmam1dt1 , Sigmam1dt2);

			T Theta = (V_dtS-V_0S)/dt; */
			
			plotfile1 << S1 << " " << S2 << " " << V_S1  << std::endl;
			//plotfile1_ref << S1 << " " << S2 << " " << V_S1_ref << " " <<  u1_ref  << " " << u2_ref << " " << u_ref  << std::endl;
            plotfile2 << S1 << " " << S2 << " " << V_S2 << " " << std::endl;//<<V_S2_ref <<  std::endl;
			plotfile3 << S1 << " " << S2 << " " << V_t <<  std::endl;
        }
        plotfile1 << std::endl;
		plotfile2 << std::endl;
		plotfile3 << std::endl;
		
    }
	time.stop();
	std::cerr << "Sensitivities computations needed " << time.elapsed() << " seconds" << std::endl;
	plotfile1.close();
	plotfile2.close();
	plotfile3.close();
	
    
}



