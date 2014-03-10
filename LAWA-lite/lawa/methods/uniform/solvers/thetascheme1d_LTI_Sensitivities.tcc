namespace lawa{

// THETASCHEME
	template<typename T, typename Basis, typename BilinearForm, typename tildeBilinearForm, typename RHSIntegral,
         typename L2ScalarProduct>
	ThetaScheme1D_LTI_Sensitivities<T, Basis, BilinearForm, tildeBilinearForm, RHSIntegral, L2ScalarProduct>::
ThetaScheme1D_LTI_Sensitivities(const T _theta, const Basis& _basis, const BilinearForm& _a,
				  const tildeBilinearForm& _tildea, RHSIntegral& _rhs,
                  const bool _time_constant_rhs,
                  const bool _use_pcg, T _assembletol, T _lintol)
    : theta(_theta), basis(_basis), tildea(_tildea),
      standardL2scalarproduct(basis), L2scalarproduct(standardL2scalarproduct),
      time_constant_rhs(_time_constant_rhs), use_pcg(_use_pcg),
      assembletol(_assembletol), lintol(_lintol),
      assembler(basis),
	  //integral(_basis, _basis),
      op_LHSMatrix(this, _a), op_RHSMatrix(this, _a), op_RHSVector(this, _rhs), prec(op_LHSMatrix),
      currentLevel(-1), P(assembler.assemblePreconditioner(prec, basis.j0))
{
	
}

	template<typename T, typename Basis, typename BilinearForm, typename tildeBilinearForm, typename RHSIntegral,
         typename L2ScalarProduct>
	ThetaScheme1D_LTI_Sensitivities<T, Basis, BilinearForm, tildeBilinearForm, RHSIntegral, L2ScalarProduct>::
ThetaScheme1D_LTI_Sensitivities(const T _theta, const Basis& _basis, const BilinearForm& _a,
				  const tildeBilinearForm& _tildea, RHSIntegral& _rhs,
                  const L2ScalarProduct& _L2scalarproduct, const bool _time_constant_rhs,
                  const bool _use_pcg, T _assembletol, T _lintol)
    : theta(_theta), basis(_basis), tildea(_tildea),
      standardL2scalarproduct(basis), L2scalarproduct(_L2scalarproduct),
      time_constant_rhs(_time_constant_rhs), use_pcg(_use_pcg),
      assembletol(_assembletol), lintol(_lintol),
      assembler(basis),
      //integral(_basis, _basis),
      op_LHSMatrix(this, _a), op_RHSMatrix(this, _a), op_RHSVector(this, _rhs), prec(op_LHSMatrix),
      currentLevel(-1), P(assembler.assemblePreconditioner(prec, basis.j0))
{
}

   
	// solve() computes a theta-step of the time-space discretization (see HilberSchwabWinter)        
template<typename T, typename Basis, typename BilinearForm, typename tildeBilinearForm, typename RHSIntegral,
         typename L2ScalarProduct>
void //flens::DenseVector<flens::Array<T> > 
ThetaScheme1D_LTI_Sensitivities<T, Basis, BilinearForm, tildeBilinearForm, RHSIntegral, L2ScalarProduct>::
solve(T time_old, T time_new, flens::DenseVector<flens::Array<T> > &u_init,
	   flens::DenseVector<flens::Array<T> > &tilde_u, int level)
{
	op_RHSVector.setTimes(time_old, time_new);
    if(level != currentLevel){
		Timer time, time2,time3;
		time.start();
		
        op_LHSMatrix.setTimes(time_old, time_new);
		
        op_RHSMatrix.setTimes(time_old, time_new);
		time3.start();
		lhsmatrix = assembler.assembleStiffnessMatrix(op_LHSMatrix, level, assembletol);
		time3.stop();
		rhsmatrix = assembler.assembleStiffnessMatrix(op_RHSMatrix, level, assembletol);
		time2.start();
		tildeA = assembler.assembleStiffnessMatrix(tildea,level, assembletol);
		time2.stop();
		P = assembler.assemblePreconditioner(prec, level);
        rhsvector = assembler.assembleRHS(op_RHSVector, level);
        currentLevel = level;
		time.stop();
		std::cout << "thetascheme1d_LTI_Sensitivities.tcc, compr-level       : " << tildea.internal_compression_level << std::endl;
		std::cout << "thetascheme1d_LTI_Sensitivities.tcc, assembling        : " << time.elapsed() << " sec" << std::endl;
		std::cout << "thetascheme1d_LTI_Sensitivities.tcc, A (lhs) assembling: " << time3.elapsed() << " sec" << std::endl;
		std::cout << "thetascheme1d_LTI_Sensitivities.tcc, tildeA assembling : " << time2.elapsed() << " sec" << std::endl;
    }
    if (!time_constant_rhs) {
        rhsvector = assembler.assembleRHS(op_RHSVector, level);
	}

 	flens::DenseVector<flens::Array<T> > u1u0 = (1-theta)* u_init;
	flens::DenseVector<flens::Array<T> > rhs = rhsmatrix * u_init + rhsvector;
	flens::DenseVector<flens::Array<T> > u(basis.mra.rangeI(level));
	
    if (use_pcg){ 
        pcg(P,lhsmatrix, u, rhs, lintol);
    }
    else {
        //int iters = gmres(lhsmatrix, u, rhs, lintol);
		int iters = pgmresm(P,lhsmatrix, u, rhs, lintol,10); 
        //std::cerr << "Solve by gmres (1) with iters = " << iters;// << std::endl;
        //pgmresm(P,lhsmatrix, u, rhs, lintol, 20);
    }
		
    u1u0 = u1u0 + theta *u;
	rhs = rhsmatrix * tilde_u - (time_new - time_old)* tildeA * u1u0;
	
	if (use_pcg){
        pcg(P,lhsmatrix, tilde_u, rhs, lintol);
    }
    else {
        int iters = pgmresm(P,lhsmatrix, tilde_u, rhs, lintol,10);
        //pgmresm(P,lhsmatrix, u, rhs, lintol, 20);
    }
	//LK -> ziemlich idiotisch, aber:
	u_init=u;
	
}

// template<typename T, typename Basis, typename BilinearForm, typename tildeBilinearForm, typename RHSIntegral,
//          typename L2ScalarProduct>
// flens::DenseVector<flens::Array<T> > 
// ThetaScheme1D_LTI_Sensitivities<T, Basis, BilinearForm, tildeBilinearForm, RHSIntegral, L2ScalarProduct>::
// solve(T time_old, T time_new, flens::DenseVector<flens::Array<T> > u_init, 
//       flens::DenseVector<flens::Array<T> > f, int level)
// {
//      if(level != currentLevel){
//          op_LHSMatrix.setTimes(time_old, time_new);
//          op_RHSMatrix.setTimes(time_old, time_new);
//          lhsmatrix = assembler.assembleStiffnessMatrix(op_LHSMatrix, level, assembletol);
//          rhsmatrix = assembler.assembleStiffnessMatrix(op_RHSMatrix, level, assembletol);
// 		 tildeA = assembler.assembleStiffnessMatrix(tildea,level, assembletol);
// 		 P = assembler.assemblePreconditioner(prec, level);
//          currentLevel = level;
//      }

	 
	 
//      flens::DenseVector<flens::Array<T> > rhs = rhsmatrix * u_init + f; //LK: aufpassen, in HilberSchwabWinter -f!
//      flens::DenseVector<flens::Array<T> > u(u_init);
//      if (use_pcg) {
//          pcg(P,lhsmatrix, u, rhs, lintol);
//      }
//      else {
//          pgmres(P,lhsmatrix, u, rhs, lintol);
//      }
//      //std::cout << cg(lhsmatrix, u, rhs) << "cg iterations" << std::endl;
//      //std::cout << pcg(P, lhsmatrix, u, rhs) << "pcg iterations" << std::endl;
//      //std::cout << "u(" << time_new << "): " << u << std::endl; 
//      return u;     
// }


template<typename T, typename Basis, typename BilinearForm, typename tildeBilinearForm, typename RHSIntegral,
         typename L2ScalarProduct>
void
ThetaScheme1D_LTI_Sensitivities<T, Basis, BilinearForm, tildeBilinearForm, RHSIntegral, L2ScalarProduct>::
setRHS(RHSIntegral& _rhs)
{
    op_RHSVector.setRHS(_rhs);
}

template<typename T, typename Basis, typename BilinearForm, typename tildeBilinearForm, typename RHSIntegral,
         typename L2ScalarProduct>
flens::SparseGeMatrix<flens::CRS<T,flens::CRS_General> > 
ThetaScheme1D_LTI_Sensitivities<T, Basis, BilinearForm, tildeBilinearForm, RHSIntegral, L2ScalarProduct>::
getLHSMatrix(int level)
{   
    if (level != currentLevel) {
        flens::SparseGeMatrix<flens::CRS<T,flens::CRS_General> > matrix
                    = assembler.assembleStiffnessMatrix(op_LHSMatrix, level);
        return matrix;
    }
    return lhsmatrix;
}

/*======================================================================================*/    
// OPERATOR_LHSMATRIX
template<typename T, typename Basis, typename BilinearForm, typename tildeBilinearForm, typename RHSIntegral,
         typename L2ScalarProduct>
ThetaScheme1D_LTI_Sensitivities<T, Basis, BilinearForm, tildeBilinearForm, RHSIntegral,L2ScalarProduct>::Operator_LHSMatrix::
Operator_LHSMatrix(ThetaScheme1D_LTI_Sensitivities<T, Basis, BilinearForm, tildeBilinearForm, RHSIntegral, L2ScalarProduct>* _scheme,
                   const BilinearForm& _a)
    : a(_a)
{   
    scheme = _scheme;
}

template<typename T, typename Basis, typename BilinearForm, typename tildeBilinearForm, typename RHSIntegral,
         typename L2ScalarProduct>
T 
ThetaScheme1D_LTI_Sensitivities<T, Basis, BilinearForm, tildeBilinearForm, RHSIntegral, L2ScalarProduct>::Operator_LHSMatrix::
operator()(XType xtype1, int j1, int k1,
           XType xtype2, int j2, int k2) const
{
    // (M + deltaT * theta * A_k+1)
    //return scheme->integral(j1, k1, xtype1, 0, j2, k2, xtype2, 0)
    //        + (time_new - time_old) * scheme->theta * a(xtype1, j1, k1, xtype2, j2, k2);
    return scheme->L2scalarproduct(xtype1, j1, k1, xtype2, j2, k2)
              + (time_new - time_old) * scheme->theta * a(xtype1, j1, k1, xtype2, j2, k2);
}


/*======================================================================================*/    
// OPERATOR_RHSMATRIX
template<typename T, typename Basis, typename BilinearForm, typename tildeBilinearForm, typename RHSIntegral,
         typename L2ScalarProduct>
ThetaScheme1D_LTI_Sensitivities<T, Basis, BilinearForm, tildeBilinearForm, RHSIntegral, L2ScalarProduct>::Operator_RHSMatrix::
Operator_RHSMatrix(const ThetaScheme1D_LTI_Sensitivities<T, Basis, BilinearForm, tildeBilinearForm, RHSIntegral,
                                           L2ScalarProduct>* _scheme,
                   const BilinearForm& _a)
    : a(_a)
{
     scheme = _scheme;
}


template<typename T, typename Basis, typename BilinearForm, typename tildeBilinearForm, typename RHSIntegral,
         typename L2ScalarProduct>
T 
ThetaScheme1D_LTI_Sensitivities<T, Basis, BilinearForm, tildeBilinearForm, RHSIntegral, L2ScalarProduct>::Operator_RHSMatrix::
operator()(XType xtype1, int j1, int k1,
           XType xtype2, int j2, int k2) const
{
   // (M - deltaT * (1-theta) * A_k)
   //return scheme->integral(j1, k1, xtype1, 0, j2, k2, xtype2, 0)
   //     - (time_new - time_old) * (1. - scheme->theta) * a(xtype1,j1,k1, xtype2, j2,k2);
    return scheme->L2scalarproduct(xtype1, j1, k1, xtype2, j2, k2)
          - (time_new - time_old) * (1. - scheme->theta) * a(xtype1,j1,k1, xtype2, j2,k2);
}

/*======================================================================================*/    
// OPERATOR_RHSVECTOR
template<typename T, typename Basis, typename BilinearForm, typename tildeBilinearForm, typename RHSIntegral,
         typename L2ScalarProduct>
ThetaScheme1D_LTI_Sensitivities<T, Basis, BilinearForm, tildeBilinearForm, RHSIntegral, L2ScalarProduct>::Operator_RHSVector::
Operator_RHSVector(const ThetaScheme1D_LTI_Sensitivities<T, Basis, BilinearForm, tildeBilinearForm, RHSIntegral,
                                           L2ScalarProduct>* _scheme,
                   RHSIntegral& _rhs)
    : rhs(_rhs)
{
     scheme = _scheme;
}

template<typename T, typename Basis, typename BilinearForm, typename tildeBilinearForm, typename RHSIntegral,
         typename L2ScalarProduct>
T 
ThetaScheme1D_LTI_Sensitivities<T, Basis, BilinearForm, tildeBilinearForm, RHSIntegral, L2ScalarProduct>::Operator_RHSVector::
operator()(XType xtype, int j, int k) const
{   
    // deltaT * (theta*f_k+1 - (1-theta)*f_k)
    return (time_new - time_old)*(scheme->theta * rhs(time_new, xtype, j, k)
                    + (1. - scheme->theta)*rhs(time_old, xtype, j, k));
} 
  
}

