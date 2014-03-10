namespace lawa {

template <typename T>
ProcessParameters2D<T,BlackScholes2D>::ProcessParameters2D
(T _r, T _sigma1, T _sigma2, T _rho, T _u11, T _u12, T _u21, T _u22)
: r(_r), sigma1(_sigma1), sigma2(_sigma2), rho(_rho), u11(_u11), u12(_u12), u21(_u21), u22(_u22)
{
	
}

template <typename T>
std::ostream& operator<<(std::ostream &s, const ProcessParameters2D<T,BlackScholes2D> &processparameters) {
    s << "_BS2D_r_" << processparameters.r      << "_sigma1_" << processparameters.sigma1
      << "_sigma2_" << processparameters.sigma2 << "_rho_" << processparameters.rho;
    return s;
}

template <typename T>
ProcessParameters2D<T,CGMYeUnivariateJump2D>::ProcessParameters2D
(T _r, T _sigma1, T _sigma2, T _rho, T _k_C1, T _k_G1, T _k_M1, T _k_Y1,
 T _k_C2, T _k_G2, T _k_M2, T _k_Y2)
: r(_r), sigma1(_sigma1), sigma2(_sigma2), rho(_rho),
  k_C1(_k_C1), k_G1(_k_G1), k_M1(_k_M1), k_Y1(_k_Y1),
  k_C2(_k_C2), k_G2(_k_G2), k_M2(_k_M2), k_Y2(_k_Y2),
  proc_param1(r,k_C1,k_G1,k_M1,k_Y1,sigma1), proc_param2(r,k_C2,k_G2,k_M2,k_Y2,sigma2)
{

}

template <typename T>
std::ostream& operator<<(std::ostream &s, const ProcessParameters2D<T,CGMYeUnivariateJump2D> &processparameters) {
    s << "_CGMYeUniv2D_r_" << processparameters.r      << "_sigma1_" << processparameters.sigma1
      << "_sigma2_" << processparameters.sigma2 << "_rho_" << processparameters.rho
      << "_C1_" << processparameters.k_C1 << "_G1_" << processparameters.k_G1
      << "_M1_" << processparameters.k_M1 << "_Y1_" << processparameters.k_Y1
      << "_C2_" << processparameters.k_C2 << "_G2_" << processparameters.k_G2
      << "_M2_" << processparameters.k_M2 << "_Y2_" << processparameters.k_Y2;
    return s;
}



	/////////////////////////////////////////////
    // //LK: Linear Transformation possible	   //
    /////////////////////////////////////////////

template <typename T>
ProcessParameters2D<T,CGMYeLinearTransformation2DJump2D>::ProcessParameters2D():
	proc_param1(0,0,0,0,0,0), proc_param2(0,0,0,0,0,0)
{
	
}
	
	
	
	
 template <typename T>
 ProcessParameters2D<T,CGMYeLinearTransformation2DJump2D>::ProcessParameters2D
 (T _r, T _sigma1, T _sigma2, T _rho, T _k_C1, T _k_G1, T _k_M1, T _k_Y1,
  T _k_C2, T _k_G2, T _k_M2, T _k_Y2, T _Sigma11, T _Sigma12, T _Sigma21, T _Sigma22 )
 : r(_r), sigma1(_sigma1), sigma2(_sigma2), rho(_rho),
   k_C1(_k_C1), k_G1(_k_G1), k_M1(_k_M1), k_Y1(_k_Y1),
   k_C2(_k_C2), k_G2(_k_G2), k_M2(_k_M2), k_Y2(_k_Y2),
   Sigma11(_Sigma11),Sigma12(_Sigma12),Sigma21(_Sigma21),Sigma22(_Sigma22), 
   proc_param1(r,k_C1,k_G1,k_M1,k_Y1,sigma1), proc_param2(r,k_C2,k_G2,k_M2,k_Y2,sigma2)
 {

 }

  template <typename T>
  std::ostream& operator<<(std::ostream &s, const ProcessParameters2D<T,CGMYeLinearTransformation2DJump2D> &processparameters) {
      s << "_CGMYeLinTrans2D_r_" << processparameters.r      << "_sigma1_" << processparameters.sigma1
        << "_sigma2_" << processparameters.sigma2 << "_rho_" << processparameters.rho
        << "_C1_" << processparameters.k_C1 << "_G1_" << processparameters.k_G1
        << "_M1_" << processparameters.k_M1 << "_Y1_" << processparameters.k_Y1
        << "_C2_" << processparameters.k_C2 << "_G2_" << processparameters.k_G2
        << "_M2_" << processparameters.k_M2 << "_Y2_" << processparameters.k_Y2
        << "_Sig11_" << processparameters.Sigma11 << "_Sig12_" << processparameters.Sigma12
        << "_Sig21_" << processparameters.Sigma21 << "_Sig22_" << processparameters.Sigma22;
      return s;
  }


template <typename T>
std::string
ProcessParameters2D<T,CGMYeLinearTransformation2DJump2D>::getSimulationName() const{
	std::stringstream s;
	s << "_CGMYeLinTrans2D_sigma1_" << sigma1
      << "_sigma2_" << sigma2 << "_rho_" << rho
      << "_C1_" << k_C1 << "_G1_" << k_G1
      << "_M1_" << k_M1 << "_Y1_" << k_Y1
      << "_C2_" << k_C2 << "_G2_" << k_G2
      << "_M2_" << k_M2 << "_Y2_" << k_Y2;
    return s.str();//s.str();	
};


	
}   // namespace lawa
