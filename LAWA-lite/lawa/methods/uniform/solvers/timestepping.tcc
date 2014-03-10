namespace lawa{
    
template <typename T, typename Solver>
TimeStepping<T,Solver>::TimeStepping(Solver& _solver, T _deltaT, int _timesteps, int _levelX)
    : solver(_solver), deltaT(_deltaT), timesteps(_timesteps), levelX(_levelX)
{    
}    
    
template <typename T, typename Solver>
flens::DenseVector<flens::Array<T> > 
TimeStepping<T,Solver>::solve(flens::DenseVector<flens::Array<T> >& u_0, bool saveSols)
{
    flens::DenseVector<flens::Array<T> > u_next(u_0.range()), u(u_0);
    if(saveSols){
        U.engine().resize(u_0.length(), timesteps+1, u_0.range().firstIndex(), 0);
        U(flens::_, 0) = u_0;
    }
    
    for(int k = 1; k <= timesteps; ++k){
        if (k%100==0) std::cerr << "TimeStepping<T,Solver>: time step = " << k << " of " << timesteps << std::endl;
        //std::cout << "k = " << k << ": u      = " << u << std::endl;
		//LK solve(T time_old, T time_new, flens::DenseVector<flens::Array<T> > u_init, 
		//LK flens::DenseVector<flens::Array<T> > f, int level)

		//LK, hier müsste f berechnet werden mit f= \tilde{A} (\theta u¹ + (1-\theta) u⁰)
		//LK, das heißt: brauchen prozedur tildeAv(v), und brauchen zweiten vector u⁰
		//LK -> neue solve- Prozedur!
        u_next = solver.solve((k-1)*deltaT, k*deltaT, u, levelX);
        //std::cout << "k = " << k << ": u_next = " << u_next << std::endl;
        u = u_next;
        if(saveSols){
            U(flens::_, k) = u;
        }
    }
    
    return u;
} 

template <typename T, typename Solver>
flens::DenseVector<flens::Array<T> > 
TimeStepping<T,Solver>::solve(flens::DenseVector<flens::Array<T> >& u_0, 
      flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >& fmatrix)
{
    assert(fmatrix.numCols() == timesteps);
    assert(fmatrix.numRows() == u_0.length());
    
    flens::DenseVector<flens::Array<T> > u_next(u_0), u(u_0);
    
    for(int k = 1; k <= timesteps; ++k){
        flens::DenseVector<flens::Array<T> > f = fmatrix(flens::_, k);
        u_next = solver.solve((k-1)*deltaT, k*deltaT, u, f, levelX);
        u = u_next;
    }
    
    return u; 
}


//solve(T time_old, T time_new, flens::DenseVector<flens::Array<T> > &u,
//	   flens::DenseVector<flens::Array<T> > &tilde_u, int level)

//LK:    this is needed to use thetatimescheme_Sensitivities
// since we need to parallel solutions, u and tilde(u)
template <typename T, typename Solver>
flens::DenseVector<flens::Array<T> > 
TimeStepping<T,Solver>::solveParallel(flens::DenseVector<flens::Array<T> >& u_0, bool saveSols)
{
    flens::DenseVector<flens::Array<T> > u_next(u_0.range()), u(u_0);
    if(saveSols){
        U.engine().resize(u_0.length(), timesteps+1, u_0.range().firstIndex(), 0);
        U(flens::_, 0) = u_0;
    }
	flens::DenseVector<flens::Array<T> > & tilde_u(u.range());
	//LK, does this set tilde_u to zero? check
	std::cout << " tilde u = " << std::cout;
	std::cout << tilde_u << std::endl;
	
    for(int k = 1; k <= timesteps; ++k){
        if (k%100==0) std::cerr << "TimeStepping<T,Solver>: time step = " << k << " of " << timesteps << std::endl;
        //std::cout << "k = " << k << ": u      = " << u << std::endl;
		//LK solve(T time_old, T time_new, flens::DenseVector<flens::Array<T> > u_init, 
		//LK flens::DenseVector<flens::Array<T> > f, int level)

		//LK, hier müsste f berechnet werden mit f= \tilde{A} (\theta u¹ + (1-\theta) u⁰)
		//LK, das heißt: brauchen prozedur tildeAv(v), und brauchen zweiten vector u⁰
		//LK -> neue solve- Prozedur!
        solver.solve((k-1)*deltaT, k*deltaT, u, tilde_u, levelX);
        //std::cout << "k = " << k << ": u_next = " << u_next << std::endl;
        
        if(saveSols){
            U(flens::_, k) = u;
        }
    }
    
    return u;
} 




template <typename T, typename Solver>
flens::DenseVector<flens::Array<T> > 
TimeStepping<T,Solver>::getResiduum(flens::DenseVector<flens::Array<T> >& u)
{
    flens::DenseVector<flens::Array<T> > Su = solve(u);
    return u - Su;
}


template <typename T, typename Solver>
flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >&
TimeStepping<T,Solver>::getSolutions()
{ 
    return U;
} 
       
template <typename T, typename Solver>
flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >&
TimeStepping<T,Solver>::getSolutions(flens::DenseVector<flens::Array<T> >& u)
{
    solve(u, true);
    return U;
}
    
} // namespace lawa

