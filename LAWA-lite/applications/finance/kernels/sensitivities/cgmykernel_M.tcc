namespace lawa {


template <typename T>
Kernel<T,CGMYeSensitivitiesM>::Kernel(const ProcessParameters1D<T,CGMY> &_params)
: params(_params), C(_params.k_C), G(_params.k_G), M(_params.k_M), Y(_params.k_Y-1)
{
	//LK Y -> Y-1 to (differentiated wrt G)
    if (((Y+1) == 0) || ((Y+1) == 1)) {
        std::cout << "The cases Y=0 and Y=1 are not implemented yet." << std::endl;
        exit(1);
    }
	
    powM_Ym5 = std::pow(M,Y-5);
    powM_Ym4 = M*powM_Ym5;
    powM_Ym3 = M*powM_Ym4;
    powM_Ym2 = M*powM_Ym3;
    powM_Ym1 = M*powM_Ym2;
    powM_Y   = M*powM_Ym1;

    powG_Ym5 = std::pow(G,Y-5);
    powG_Ym4 = G*powG_Ym5;
    powG_Ym3 = G*powG_Ym4;
    powG_Ym2 = G*powG_Ym3;
    powG_Ym1 = G*powG_Ym2;
    powG_Y   = G*powG_Ym1;

    CdivY    = C/Y;
    Cdiv1mY  = C/(1.-Y);
    Mdiv1mY  = M/(1.-Y);
    Gdiv1mY  = G/(1.-Y);
    OnedivSix= 1./6.;

	//LK, M-sensitivity limits from above are zero
    constants[0] = - 0.5L     * C * powM_Ym2 * boost::math::tgamma(2-Y); // c_3(O+) was +
    constants[1] =  0.0; // -0.5L     * C * powG_Ym2 * boost::math::tgamma(2-Y); // c_3(0-)
    constants[2] = - 1.L/6.L  * C * powM_Ym3 * boost::math::tgamma(3-Y); // c_4(0+) was +
    constants[3] =  0.0; // LK, was 1.L/6.L  * C * powG_Ym3 * boost::math::tgamma(3-Y); // c_4(0-)
    constants[4] = - 1.L/24.L * C * powM_Ym4 * boost::math::tgamma(4-Y); // c_5(0+) was +
    constants[5] =  0.0; //LK, was -1.L/24.L * C * powG_Ym4 * boost::math::tgamma(4-Y); // c_5(0-)
    constants[6] = - 1.L/120.L* C * powM_Ym5 * boost::math::tgamma(5-Y); // c_6(0+) was +
    constants[7] =  0.0;//LK, was 1.L/120.L* C * powG_Ym5 * boost::math::tgamma(5-Y); // c_6(0-)
    //LK, was constants[8] =  nthTailIntegral(10e-16,7,ZeroAtInfinity);
    //LK, was constants[9] =  nthTailIntegral(-10e-16,7,ZeroAtInfinity);

	//LK for G-derivative -> use new \tilde{k}^{-n} := 1_{x<0} (k^{-n}(x) - k^n(0-))
	//                       and c_n := - k^{-n}(0-)
    
	//LK new constants c_n := - k^{-n}(0-) = \tilde{k}^{-n}(0+)
    //LK, was c3 = constants[0]-constants[1];
	
	c3 = constants[0]-constants[1];
    //LK, was c4 = constants[2]-constants[3];
	c4 =  constants[2]- constants[3];
    //LK, was c5 = constants[4]-constants[5];
	c5 = constants[4]- constants[5];
	//LK, was c6 = constants[6]-constants[7];
	c6 = constants[6] -constants[7];

    // c3 = constants[0]-constants[1];
    // c4 = constants[2]-constants[3];
    // c5 = constants[4]-constants[5];
    // c6 = constants[6]-constants[7];

    ExpXmOnemX_k_pos = C*boost::math::tgamma(-Y)*( std::pow(M-1,Y)-std::pow(M,Y)+Y*std::pow(M,Y-1));
    ExpXmOnemX_k_neg = C*boost::math::tgamma(-Y)*( std::pow(G+1,Y)-std::pow(G,Y)-Y*std::pow(G,Y-1));
    ExpXmOnemX_k     = ExpXmOnemX_k_pos + ExpXmOnemX_k_neg;
    ExpXmOne_k1_pos =  ExpXmOnemX_k_pos;
    ExpXmOne_k1_neg =  ExpXmOnemX_k_neg;
    ExpXmOne_k2_pos =  (CdivY/(1-Y))*(boost::math::tgamma(2-Y)*( powM_Y-2*powM_Ym1-std::pow(M-1,Y) )
                                      + boost::math::tgamma(3-Y)*powM_Ym1  ) - constants[0];
    ExpXmOne_k2_neg = -(CdivY/(1-Y))*(boost::math::tgamma(2-Y)*(-powG_Y-2*powG_Ym1+std::pow(G+1,Y) )
                                      + boost::math::tgamma(3-Y)*powG_Ym1  ) + constants[1];


/*
    std::cout << "CGMY-Utils: c3=" << c3 << ", c4=" << c4 << ", c5=" << c5 << ", c6=" << c6 << std::endl;
    std::cout << "CGMY-Utils: constants[0]=" << constants[0] << ", constants[1]=" << constants[1] << ", constants[2]=" << constants[2] << ", constants[3]=" << constants[3] << std::endl;
    std::cout << "CGMY-Utils: constants[4]=" << constants[4] << ", constants[5]=" << constants[5] << ", constants[6]=" << constants[6] << ", constants[7]=" << constants[7] << std::endl;

    std::cout << "CGMY-Utils: ExpXmOnemX_k_pos = " << ExpXmOnemX_k_pos << ", ExpXmOnemX_k_neg = " << ExpXmOnemX_k_neg << std::endl;
    std::cout << "CGMY-Utils: ExpXmOne_k1_pos = " << ExpXmOne_k1_pos << ", ExpXmOne_k1_neg = " << ExpXmOne_k1_neg << std::endl;
    std::cout << "CGMY-Utils: ExpXmOne_k2_pos = " << ExpXmOne_k2_pos << ", ExpXmOne_k2_neg = " << ExpXmOne_k2_neg << std::endl;
    std::cout << "CGMY-Utils: Initialization finished." << std::endl;
*/

}

template <typename T>
Kernel<T,CGMYeSensitivitiesM>::Kernel(const ProcessParameters1D<T,CGMYe> &_params)
  : params(_params.r, _params.k_C, _params.k_G, _params.k_M, _params.k_Y-1),//-1
	C(_params.k_C), G(_params.k_G), M(_params.k_M), Y(_params.k_Y-1)//-1
{
	//LK Y -> Y-1
    if (((Y+1) == 0) || ((Y+1) == 1)) {
        std::cout << "The cases Y=0 and Y=1 are not implemented yet." << std::endl;
        exit(1);
    }

    powM_Ym5 = std::pow(M,Y-5);
    powM_Ym4 = M*powM_Ym5;
    powM_Ym3 = M*powM_Ym4;
    powM_Ym2 = M*powM_Ym3;
    powM_Ym1 = M*powM_Ym2;
    powM_Y   = M*powM_Ym1;

    powG_Ym5 = std::pow(G,Y-5);
    powG_Ym4 = G*powG_Ym5;
    powG_Ym3 = G*powG_Ym4;
    powG_Ym2 = G*powG_Ym3;
    powG_Ym1 = G*powG_Ym2;
    powG_Y   = G*powG_Ym1;

    CdivY    = C/Y;
    Cdiv1mY  = C/(1.L-Y);
    Mdiv1mY  = M/(1.L-Y);
    Gdiv1mY  = G/(1.L-Y);
    OnedivSix= 1.L/6.L;

	/*
	//LK, G-sensitivity limits from above are zero
    constants[0] =  0.0; //LK, was 0.5L     * C * powM_Ym2 * boost::math::tgamma(2-Y);
    constants[1] = 0.5L     * C * powG_Ym2 * boost::math::tgamma(2-Y); //LK, was -
    constants[2] =  0.0; //LK, was 1.L/6.L  * C * powM_Ym3 * boost::math::tgamma(3-Y);
    constants[3] = - 1.L/6.L  * C * powG_Ym3 * boost::math::tgamma(3-Y); //LK, was +
    constants[4] =  0.0; //LK, was 1.L/24.L * C * powM_Ym4 * boost::math::tgamma(4-Y);
    constants[5] = 1.L/24.L * C * powG_Ym4 * boost::math::tgamma(4-Y); //LK, was -
    constants[6] =  0.0; //LK, was 1.L/120.L* C * powM_Ym5 * boost::math::tgamma(5-Y);
    constants[7] = - 1.L/120.L* C * powG_Ym5 * boost::math::tgamma(5-Y); //LK, was +
    //LK, was constants[8] =  nthTailIntegral(10e-16,7,ZeroAtInfinity);
    //LK, was constants[9] =  nthTailIntegral(-10e-16,7,ZeroAtInfinity);

	//LK for G-derivative -> use new \tilde{k}^{-n} := 1_{x<0} (k^{-n}(x) - k^n(0-))
	//                       and c_n := - k^{-n}(0-)

	//LK new constants c_n := - k^{-n}(0-)
    //LK, was c3 = constants[0]-constants[1];
	// (-1)*(-1)
	c3 = -constants[1];
    //LK, was c4 = constants[2]-constants[3];
	c4 = -constants[3];
    //LK, was c5 = constants[4]-constants[5];
	c5 = -constants[5];
	//LK, was c6 = constants[6]-constants[7];
	c6 = -constants[7];
	*/

		//LK, M-sensitivity limits from above are zero
    constants[0] = - 0.5L     * C * powM_Ym2 * boost::math::tgamma(2-Y); // c_3(O+) was +
    constants[1] =  0.0; // -0.5L     * C * powG_Ym2 * boost::math::tgamma(2-Y); // c_3(0-)
    constants[2] = - 1.L/6.L  * C * powM_Ym3 * boost::math::tgamma(3-Y); // c_4(0+) was +
    constants[3] =  0.0; // LK, was 1.L/6.L  * C * powG_Ym3 * boost::math::tgamma(3-Y); // c_4(0-)
    constants[4] = - 1.L/24.L * C * powM_Ym4 * boost::math::tgamma(4-Y); // c_5(0+) was +
    constants[5] =  0.0; //LK, was -1.L/24.L * C * powG_Ym4 * boost::math::tgamma(4-Y); // c_5(0-)
    constants[6] = - 1.L/120.L* C * powM_Ym5 * boost::math::tgamma(5-Y); // c_6(0+) was +
    constants[7] =  0.0;//LK, was 1.L/120.L* C * powG_Ym5 * boost::math::tgamma(5-Y); // c_6(0-)
    //LK, was constants[8] =  nthTailIntegral(10e-16,7,ZeroAtInfinity);
    //LK, was constants[9] =  nthTailIntegral(-10e-16,7,ZeroAtInfinity);

	//LK for G-derivative -> use new \tilde{k}^{-n} := 1_{x<0} (k^{-n}(x) - k^n(0-))
	//                       and c_n := - k^{-n}(0-)
    
	//LK new constants c_n := - k^{-n}(0-) = \tilde{k}^{-n}(0+)
    //LK, was c3 = constants[0]-constants[1];
	
	c3 = constants[0]-constants[1];
    //LK, was c4 = constants[2]-constants[3];
	c4 =  constants[2]- constants[3];
    //LK, was c5 = constants[4]-constants[5];
	c5 = constants[4]- constants[5];
	//LK, was c6 = constants[6]-constants[7];
	c6 = constants[6] -constants[7];
	

	//LK like in Tops: page 76 
	ExpXmOnemX_k_pos = C*boost::math::tgamma(-Y)*( std::pow(M-1,Y)-std::pow(M,Y)+Y*std::pow(M,Y-1));
    ExpXmOnemX_k_neg = C*boost::math::tgamma(-Y)*( std::pow(G+1,Y)-std::pow(G,Y)-Y*std::pow(G,Y-1));
    ExpXmOnemX_k     = ExpXmOnemX_k_pos + ExpXmOnemX_k_neg;
    ExpXmOne_k1_pos =  ExpXmOnemX_k_pos;
    ExpXmOne_k1_neg =  ExpXmOnemX_k_neg;
    ExpXmOne_k2_pos =  (CdivY/(1-Y))*(boost::math::tgamma(2-Y)*( powM_Y-2*powM_Ym1-std::pow(M-1,Y) )
                                      + boost::math::tgamma(3-Y)*powM_Ym1  ) - constants[0];
    ExpXmOne_k2_neg = -(CdivY/(1-Y))*(boost::math::tgamma(2-Y)*(-powG_Y-2*powG_Ym1+std::pow(G+1,Y) )
                                      + boost::math::tgamma(3-Y)*powG_Ym1  ) + constants[1];


/*
    std::cout << "CGMY-Utils: c3=" << c3 << ", c4=" << c4 << ", c5=" << c5 << ", c6=" << c6 << std::endl;
    std::cout << "CGMY-Utils: constants[0]=" << constants[0] << ", constants[1]=" << constants[1] << ", constants[2]=" << constants[2] << ", constants[3]=" << constants[3] << std::endl;
    std::cout << "CGMY-Utils: constants[4]=" << constants[4] << ", constants[5]=" << constants[5] << ", constants[6]=" << constants[6] << ", constants[7]=" << constants[7] << std::endl;

    std::cout << "CGMY-Utils: ExpXmOnemX_k_pos = " << ExpXmOnemX_k_pos << ", ExpXmOnemX_k_neg = " << ExpXmOnemX_k_neg << std::endl;
    std::cout << "CGMY-Utils: ExpXmOne_k1_pos = " << ExpXmOne_k1_pos << ", ExpXmOne_k1_neg = " << ExpXmOne_k1_neg << std::endl;
    std::cout << "CGMY-Utils: ExpXmOne_k2_pos = " << ExpXmOne_k2_pos << ", ExpXmOne_k2_neg = " << ExpXmOne_k2_neg << std::endl;
    std::cout << "CGMY-Utils: Initialization finished." << std::endl;
*/

}


	//LK times (-1), 1_{x<0}
template <typename T>
T
Kernel<T,CGMYeSensitivitiesM>::TailIntegralMoments(T x, int l) const {
    if (l == 0) {
        assert(fabs(x)>0);
        if (x>0) return -(C/Y)*( (std::pow(x,-Y)*exp(-M*x)) * (1+(M/(1-Y))*x)
								- std::pow(M,Y)/(1-Y)*boost::math::tgamma(2-Y,M*x)); //LK was +
        else     return  0.0; /*-(C/Y)*( (std::pow(-x,-Y)*exp(G*x)) * (1+(G/(1-Y))*(-x))
						   - std::pow(G,Y)/(1-Y)*boost::math::tgamma(2-Y,-G*x)); // was +*/
    }
    else if (l == 1) {
        if (x>0) return - C/(1-Y)*( -std::pow(x,1-Y)*exp(-M*x)
								   +std::pow(M,Y-1)*boost::math::tgamma(2-Y,M*x) ); //LK, was +
        else     return  0.0;/*- C/(1-Y)*( -std::pow(-x,1-Y)*exp(G*x)
						   +std::pow(G,Y-1)*boost::math::tgamma(2-Y,-G*x) ); // was + */
    }
    else if (l >= 2) {
        if (x>0) return - C*std::pow(M,Y-l)*boost::math::tgamma(l-Y,M*x); //LK, was +
        else     return  0.0;/*-  C*std::pow(G,Y-l)*boost::math::tgamma(l-Y,-G*x); // was +*/
    }
    else exit(1);
}

//LK, not needed for sensitivities wrt G
// template <typename T>
// T
// Kernel_G<T,CGMY>::nthTailIntegral(T x, int n, AntiDerivativeType type) const {
//     if (n == 0 && type!=ZeroAtInfinity) {
//         std::cout << "n=0: type should be ZeroAtInfinity." << std::endl;
//         exit(1);
//     }
//     if (n == 1 && Y > 1 && type!=ZeroAtInfinity) {
//         std::cout << "n=1, Y>1: type should be ZeroAtInfinity." << std::endl;
//         exit(1);
//     }
//     if (type == ZeroAtZero && n>7) {
//         std::cout << "For type=ZeroAtZero, n should be <= 7." << std::endl;
//         exit(1);
//     }

//     T sum = 0.0;
//     if (type == ZeroAtZero && x == 0) return 0.0;
//     //if (type == ZeroAtInfinity && n == 2 && fabs(x)<1e-13) return 0.;
//     if (x > 0) {
//         for (int k=0; k<=n-1; ++k) {
//             if ((n-k-1) % 2 == 0) sum += boost::math::binomial_coefficient<T>(n-1,k)
//                                          *std::pow(x,n-1-k)*TailIntegralMoments(x,k);
//             else                  sum -= boost::math::binomial_coefficient<T>(n-1,k)
//                                          *std::pow(x,n-1-k)*TailIntegralMoments(x,k);
//         }
//         if (type==ZeroAtInfinity) return (1.0/boost::math::factorial<T>(n-1))*sum;
//         else                      return (1.0/boost::math::factorial<T>(n-1))*sum
//                                          - constants[2*(n-3)];
//     }
//     else {
//         for (int k=0; k<=n-1; ++k) {
//             if ((n-k-1) % 2 == 0)   sum += boost::math::binomial_coefficient<T>(n-1,k)
//                                            *std::pow(-x,n-1-k)*TailIntegralMoments(x,k);
//             else                    sum -= boost::math::binomial_coefficient<T>(n-1,k)
//                                            *std::pow(-x,n-1-k)*TailIntegralMoments(x,k);
//         }
//         if (type==ZeroAtInfinity) {
//             if (n % 2 == 0) return  (1.0/boost::math::factorial<T>(n-1))*sum;
//             else            return -(1.0/boost::math::factorial<T>(n-1))*sum;
//         }
//         else {
//             if (n % 2 == 0) return  (1.0/boost::math::factorial<T>(n-1))*sum - constants[2*(n-3)+1];
//             else            return -(1.0/boost::math::factorial<T>(n-1))*sum - constants[2*(n-3)+1];
//         }
//     }
// }


//LK, changed
template <typename T>
T
Kernel<T,CGMYeSensitivitiesM>::FirstTailIntegral(T x) const {
    //return this->nthTailIntegral(x, 1, ZeroAtInfinity);

    if (x>0)    return - (C/Y)*( (std::pow(x,-Y)*exp(-M*x)) * (1+(M/(1-Y))*x)
								 - std::pow(M,Y)/(1-Y)*boost::math::tgamma(2-Y,M*x)); //LK, was +
    else        return 0.0;/*(C/Y)*( (std::pow(-x,-Y)*exp(G*x)) * (1+(G/(1-Y))*(-x))
							   - std::pow(G,Y)/(1-Y)*boost::math::tgamma(2-Y,-G*x)); //LK was - */
}


// template <typename T>
// T
// Kernel_G<T,CGMY>::FirstTailFirstExpMomIntegral(T x) const {
//     if (x>0)    return  (C/Y)*( (std::pow(x,-Y)*exp(-(M-1)*x)) * (1+((M-1)/(1-Y))*x)
//                                 - std::pow(M-1,Y)/(1-Y)*boost::math::tgamma(2-Y,(M-1)*x)
//                               );
//     else        return -(C/Y)*( (std::pow(-x,-Y)*exp((G+1)*x)) * (1+((G+1)/(1-Y))*(-x))
//                                 - std::pow(G+1,Y)/(1-Y)*boost::math::tgamma(2-Y,-(G+1)*x)
//                               );
// }


	//LK, changed
template <typename T>
T
Kernel<T,CGMYeSensitivitiesM>::SecondTailIntegral(T x) const{
    //return this->nthTailIntegral(x, 2, ZeroAtInfinity);
    assert(fabs(x)>0);
	//LK, tmp!!!!!	
    if (x>0) return  -(C/(Y*(1-Y)))*(  boost::math::tgamma(2-Y,M*x)*(2.0L*std::pow(M,Y-1)+std::pow(M,Y)*x )
									   	-boost::math::tgamma(3-Y,M*x)*std::pow(M,Y-1)
									  -std::pow(x,1-Y)*exp(-M*x) ); //LK, was +
    else  {
        x=fabs(x);
		return   0.0;/*-(C/(Y*(1-Y)))*(  boost::math::tgamma(2-Y,G*x)*(
                         2.0L*std::pow(G,Y-1)+std::pow(G,Y)*x )
                        -boost::math::tgamma(3-Y,G*x)*std::pow(G,Y-1)
								   -std::pow(x,1-Y)*exp(-G*x) ); //LK, was + */
    }
}

// emplate <typename T>
// T
// Kernel_G<T,CGMY>::SecondTailFirstExpMomIntegral(T x) const{
//     return FirstTailFirstExpMomIntegral(x) - exp(x)*FirstTailIntegral(x);
// }


//LK, changed
 template <typename T>
 T
 Kernel<T,CGMYeSensitivitiesM>::ThirdTailIntegral(T x) const {
     //return this->nthTailIntegral(x, 3, ZeroAtInfinity);
     assert(fabs(x)>0);
     if (x>0) return - C/(Y*(1-Y))*( - boost::math::tgamma(2-Y,M*x)*( 2*std::pow(M,Y-1)*x
																	 +0.5L*std::pow(M,Y)*x*x+std::pow(M,Y-2) )
									+ boost::math::tgamma(3-Y,M*x)*( 2*std::pow(M,Y-2)+x*std::pow(M,Y-1) )
									+ boost::math::tgamma(4-Y,M*x)*( -0.5L*std::pow(M,Y-2) )   ); //LK, was +
     else {
         x = fabs(x);
		 return 0.0; /*C/(Y*(1-Y))*( - boost::math::tgamma(2-Y,G*x)*( 2*std::pow(G,Y-1)*x
																+0.5*std::pow(G,Y)*x*x+std::pow(G,Y-2) )
							   + boost::math::tgamma(3-Y,G*x)*( 2*std::pow(G,Y-2)+x*std::pow(G,Y-1) )
							   + boost::math::tgamma(4-Y,G*x)*( -0.5*std::pow(G,Y-2) )   ); //LK, was - */
     }
 }

// template <typename T>
// T
// Kernel_G<T,CGMY>::ThirdTailFirstExpMomIntegral(T x) const{
//     //return SecondTailFirstExpMomIntegral(x) - exp(x)* SecondTailIntegral(x);
//     if (x>0) {
//         T MX    = M*x;
//         T ExpMMX = std::exp(-MX);
//         T PowMX_2mY = std::pow(MX,2-Y);
//         T tgamma_2mY_MX = boost::math::tgamma(2-Y,MX);
//         T tgamma_3mY_MX = PowMX_2mY*ExpMMX + tgamma_2mY_MX*(2-Y);
//         T SecondTailIntegralPart           =  (C/(Y*(1-Y)))*(
//                                                     tgamma_2mY_MX*( 2.0*powM_Ym1 + powM_Y*x )
//                                                    -tgamma_3mY_MX*powM_Ym1-std::pow(x,1-Y)*ExpMMX
//                                                              );
//         T FirstTailFirstExpMomIntegralPart = (C/Y)*(  (std::pow(x,-Y)*exp(-(M-1)*x))
//                                                     * (1+((M-1)/(1-Y))*x)
//                                            -std::pow(M-1,Y)/(1-Y)*boost::math::tgamma(2-Y,(M-1)*x));
//         T FirstTailIntegralPart            = (C/Y)*( (std::pow(x,-Y)*ExpMMX) * (1+(M/(1-Y))*x)
//                                              - powM_Y/(1-Y)*tgamma_2mY_MX) ;
//         T ret = -exp(x)*SecondTailIntegralPart
//                 + (FirstTailFirstExpMomIntegralPart - exp(x)*FirstTailIntegralPart);

//         return ret;
//     }
//     else  {
//         T X = -x;
//         T GX    = G*X;
//         T ExpGMX = std::exp(-GX);
//         T PowGX_2mY = std::pow(GX,2-Y);
//         T tgamma_2mY_GX = boost::math::tgamma(2-Y,GX);
//         T tgamma_3mY_GX = PowGX_2mY*ExpGMX + tgamma_2mY_GX*(2-Y);
//         T SecondTailIntegralPart            = (C/(Y*(1-Y)))*(
//                                                 tgamma_2mY_GX*( 2.0*powG_Ym1 + powG_Y*X )
//                                                -tgamma_3mY_GX*powG_Ym1 -std::pow(X,1-Y)*ExpGMX
//                                                              );
//         T FirstTailFirstExpMomIntegralPart =-(C/Y)*( (std::pow(-x,-Y)*exp((G+1)*x))
//                                                     *(1+((G+1)/(1-Y))*(-x))
//                                                     -std::pow(G+1,Y)/(1-Y)
//                                                     *boost::math::tgamma(2-Y,-(G+1)*x)
//                                                     );
//         T FirstTailIntegralPart             =-(C/Y)*( (std::pow(-x,-Y)*ExpGMX)
//                                                      *(1+(G/(1-Y))*(-x))
//                                                      -powG_Y/(1-Y)*tgamma_2mY_GX);

//         T ret = -exp(x)*SecondTailIntegralPart + (FirstTailFirstExpMomIntegralPart
//                 -exp(x)*FirstTailIntegralPart);
//         return ret;
//     }
// }


//LK, changed
template <typename T>
T
Kernel<T,CGMYeSensitivitiesM>::ForthTailIntegral(T x) const {
    //return this->nthTailIntegral(x, 4, ZeroAtInfinity);
    assert(fabs(x)>0);

    if (x>0) {
        T MX    = M*x;
        T ExpMMX = std::exp(-MX);
        T PowMX_2mY = std::pow(MX,2-Y);
        T tgamma_2mY_MX = boost::math::tgamma(2-Y,MX);
        T tgamma_3mY_MX = PowMX_2mY*ExpMMX + tgamma_2mY_MX*(2-Y);
        T PowX_3mY = std::pow(x,3-Y);

        return -( -PowX_3mY*ExpMMX*( OnedivSix*CdivY*(1 + Mdiv1mY*x) + 0.5*Cdiv1mY )
				 +OnedivSix * C * ( powM_Ym3*tgamma_3mY_MX )
				 -0.5 * C * ( x*powM_Ym2*tgamma_2mY_MX)
				 +0.5 * Cdiv1mY * ( x*x*powM_Ym1*tgamma_2mY_MX  )
				  +OnedivSix * CdivY * ( x*x*x*powM_Y*tgamma_2mY_MX/(1-Y))); //LK, was +
    }
    else  {/*
        T X     = -x;
        T GX    = G*X;
        T ExpGMX = std::exp(-GX);
        T PowGX_2mY = std::pow(GX,2-Y);
        T tgamma_2mY_GX = boost::math::tgamma(2-Y,GX);
        T tgamma_3mY_GX = PowGX_2mY*ExpGMX + tgamma_2mY_GX*(2-Y);
        T PowX_3mY = std::pow(X,3-Y);*/

        return 0.0; /*-(- PowX_3mY*ExpGMX*( OnedivSix*CdivY*(1 + Gdiv1mY*X) + 0.5*Cdiv1mY )
				 +OnedivSix * C * ( powG_Ym3*tgamma_3mY_GX )
				 -0.5 * C * ( X*powG_Ym2*tgamma_2mY_GX)
				 +0.5 * Cdiv1mY * ( X*X*powG_Ym1*tgamma_2mY_GX  )
				 +OnedivSix * CdivY * ( X*X*X*powG_Y*tgamma_2mY_GX/(1-Y)) );//LK, was + */
    }
}

// template <typename T>
// T
// Kernel_G<T,CGMY>::ForthTailFirstExpMomIntegral(T x) const{
//     return ThirdTailFirstExpMomIntegral(x) - exp(x)*ThirdTailIntegral(x);
// }


//LK, changed
template <typename T>
T
Kernel<T,CGMYeSensitivitiesM>::FifthTailIntegral(T x) const {
	//return this->nthTailIntegral(x, 5, ZeroAtInfinity);
	assert(fabs(x)>0);

	if (x>0) {
		T MX    = M*x;
		T ExpMMX = std::exp(-MX);
		T PowMX_3mY = std::pow(MX,3-Y);
		T PowMX_2mY = PowMX_3mY / MX;
		T tgamma_2mY_MX = boost::math::tgamma(2-Y,MX);
		T tgamma_3mY_MX = PowMX_2mY*ExpMMX + tgamma_2mY_MX*(2-Y);
		T tgamma_4mY_MX = PowMX_3mY*ExpMMX + tgamma_3mY_MX*(3-Y);
		T PowX_4mY = std::pow(x,4-Y);

		return -(PowX_4mY*ExpMMX*( CdivY*(1 + Mdiv1mY*x)/24 + Cdiv1mY/6 )
				 -CdivY/24 *     ( x*x*x*x*powM_Y*tgamma_2mY_MX/(1-Y))
				 -Cdiv1mY/6 *    ( x*x*x*powM_Ym1*tgamma_2mY_MX  )
				 +C/4 *          ( x*x*powM_Ym2*tgamma_2mY_MX)
				 -C/6 *          ( x*powM_Ym3*tgamma_3mY_MX )
				 +C/24*          ( powM_Ym4*tgamma_4mY_MX)); // LK, was +
	}
	else  {/*
		T X = -x;
		T GX    = G*X;
		T ExpGMX = std::exp(-GX);
		T PowGX_3mY = std::pow(GX,3-Y);
		T PowGX_2mY = PowGX_3mY / GX;
		T tgamma_2mY_GX = boost::math::tgamma(2-Y,GX);
		T tgamma_3mY_GX = PowGX_2mY*ExpGMX + tgamma_2mY_GX*(2-Y);
		T tgamma_4mY_GX = PowGX_3mY*ExpGMX + tgamma_3mY_GX*(3-Y);
		T PowX_4mY = std::pow(X,4-Y); */

		return 0.0; /*(PowX_4mY*ExpGMX*( CdivY*(1 + Gdiv1mY*X)/24 + Cdiv1mY/6 )
				-CdivY/24 *   ( X*X*X*X*powG_Y*tgamma_2mY_GX/(1-Y))
				-Cdiv1mY/6 *  ( X*X*X*powG_Ym1*tgamma_2mY_GX  )
				+C/4 *        ( X*X*powG_Ym2*tgamma_2mY_GX)
				-C/6 *        ( X*powG_Ym3*tgamma_3mY_GX )
				+C/24*        ( powG_Ym4*tgamma_4mY_GX)); //LK, was */
	}
}


//LK, changed
template <typename T>
T
Kernel<T,CGMYeSensitivitiesM>::SixthTailIntegral(T x) const {
    //return this->nthTailIntegral(x, 6, ZeroAtInfinity);
    assert(fabs(x)>0);

    if (x>0) {
        T MX    = M*x;
        T ExpMMX = std::exp(-MX);
        T PowMX_4mY = std::pow(MX,4-Y);
        T PowMX_3mY = PowMX_4mY / MX;
        T PowMX_2mY = PowMX_3mY / MX;
        T tgamma_2mY_MX = boost::math::tgamma(2-Y,MX);
        T tgamma_3mY_MX = PowMX_2mY*ExpMMX + tgamma_2mY_MX*(2-Y);
        T tgamma_4mY_MX = PowMX_3mY*ExpMMX + tgamma_3mY_MX*(3-Y);
        T tgamma_5mY_MX = PowMX_4mY*ExpMMX + tgamma_4mY_MX*(4-Y);

        T PowX_5mY = std::pow(x,5-Y);

		return -(-PowX_5mY*ExpMMX*( CdivY*(1 + Mdiv1mY*x)/120 + Cdiv1mY/24 )
                 +C/120*        ( powM_Ym5*tgamma_5mY_MX) +
                 x*             ( -C/24 *powM_Ym4*tgamma_4mY_MX +
								  x*( C/12*powM_Ym3*tgamma_3mY_MX +
									  x*(-C/12*powM_Ym2*tgamma_2mY_MX +
										 x*( Cdiv1mY/24*powM_Ym1*tgamma_2mY_MX +
											 x*( CdivY/120*powM_Y*tgamma_2mY_MX/(1-Y) ) ) ) ) ) ); //LK, was +
    }
    else  {/*
        T X = -x;
        T GX    = G*X;
        T ExpGMX = std::exp(-GX);
        T PowGX_4mY = std::pow(GX,4-Y);
        T PowGX_3mY = PowGX_4mY / GX;
        T PowGX_2mY = PowGX_3mY / GX;
        T tgamma_2mY_GX = boost::math::tgamma(2-Y,GX);
        T tgamma_3mY_GX = PowGX_2mY*ExpGMX + tgamma_2mY_GX*(2-Y);
        T tgamma_4mY_GX = PowGX_3mY*ExpGMX + tgamma_3mY_GX*(3-Y);
        T tgamma_5mY_GX = PowGX_4mY*ExpGMX + tgamma_4mY_GX*(4-Y);
        T PowX_5mY = std::pow(X,5-Y);

        return (PowX_5mY*ExpGMX*( CdivY*(1 + Gdiv1mY*X)/120 + Cdiv1mY/24 )
                +C/120*        ( powG_Ym5*tgamma_5mY_GX) +
                X*             ( -C/24 *powG_Ym4*tgamma_4mY_GX +
                X*             ( C/12*powG_Ym3*tgamma_3mY_GX +
                X*             (-C/12*powG_Ym2*tgamma_2mY_GX +
                X*             ( Cdiv1mY/24*powG_Ym1*tgamma_2mY_GX +
				X*             ( CdivY/120*powG_Y*tgamma_2mY_GX/(1-Y) ) ) ) ) ) ); //LK, was -*/
		return 0.0;
    }
}



}   // namespace lawa
