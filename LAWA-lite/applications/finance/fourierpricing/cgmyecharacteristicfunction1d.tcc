namespace lawa {

template <typename T>
CharacteristicFunction1D<T,CGMYe>::CharacteristicFunction1D
                                   (const ProcessParameters1D<T,CGMYe> &_processparameters)
 : processparameters(_processparameters), drift(0.)
{
    T C     = processparameters.k_C;
    T G     = processparameters.k_G;
    T M     = processparameters.k_M;
    T Y     = processparameters.k_Y;
    T sigma = processparameters.sigma;
    drift = -C*boost::math::tgamma(-Y)*(std::pow(M-1,Y) - std::pow(M,Y) + std::pow(G+1,Y)
            - std::pow(G,Y))  - 0.5*sigma*sigma;;
	//std::cout << "drift (cgmyechar1) = " << drift << std::endl; ;
	//std::getchar();
}

// 	//LK: neu
// 	template <typename T>
// CharacteristicFunction1D<T,CGMYe>::CharacteristicFunction1D
// 	(T _r,T _C1,T _C2, T _C, T _G, T _M, T _Y, T _sigma1, T _sigma2, T _sigma, T _drift)
// 		: processparameters(_r,_C,_G,_M,_Y,_sigma), drift(_drift)
// {
// 	    T drift1 = -_C1*boost::math::tgamma(-_Y)*(std::pow(_M-1,_Y) - std::pow(_M,_Y) + std::pow(_G+1,_Y)
//             - std::pow(_G,_Y))  - 0.5*_sigma1*_sigma1;
// 		T drift2 = -_C2*boost::math::tgamma(-_Y)*(std::pow(_M-1,_Y) - std::pow(_M,_Y) + std::pow(_G+1,_Y)
//             - std::pow(_G,_Y))  - 0.5*_sigma2*_sigma2;
// 		drift = _drift; //drift1 + drift2 +_r;
// 		std::cout << "drift (cgmyechar2) = " << drift <<
// 			"drift 1 = " << drift1 << "drift 2 = " << drift2 <<std::endl;;
//     // T C     = processparameters.k_C;
//     // T G     = processparameters.k_G;
//     // T M     = processparameters.k_M;
//     // T Y     = processparameters.k_Y;
//     // T sigma = processparameters.sigma;
//     // drift = -C*boost::math::tgamma(-Y)*(std::pow(M-1,Y) - std::pow(M,Y) + std::pow(G+1,Y)
//     //         - std::pow(G,Y))  - 0.5*sigma*sigma;;
// }

	
template <typename T>
gsl_complex
CharacteristicFunction1D<T,CGMYe>::operator()(T t, T v)
{
    T C     = processparameters.k_C;
    T G     = processparameters.k_G;
    T M     = processparameters.k_M;
    T Y     = processparameters.k_Y;
    T sigma = processparameters.sigma;

    gsl_complex z1 = gsl_complex_rect(M,-v);
    gsl_complex z2 = gsl_complex_rect(G, v);
    gsl_complex z1_Y = gsl_complex_pow_real(z1, Y);
    gsl_complex z2_Y = gsl_complex_pow_real(z2, Y);

    gsl_complex tmp = gsl_complex_add(z1_Y,z2_Y);
    tmp = gsl_complex_add_real(tmp,-std::pow(M,Y)-std::pow(G,Y));
    tmp = gsl_complex_mul_real(tmp, C*t*boost::math::tgamma(-Y));
    tmp = gsl_complex_add(tmp, gsl_complex_rect(0.0,v*t*drift));
    tmp = gsl_complex_add_real(tmp,-0.5*sigma*sigma*v*v);
    tmp = gsl_complex_exp(tmp);

    return tmp;
}

template <typename T>
gsl_complex
CharacteristicFunction1D<T,CGMYe>::operator()(T t, gsl_complex v)
{
    T C     = processparameters.k_C;
    T G     = processparameters.k_G;
    T M     = processparameters.k_M;
    T Y     = processparameters.k_Y;
    T sigma = processparameters.sigma;

    gsl_complex z1 = gsl_complex_rect(M+v.dat[1],-v.dat[0]);
    gsl_complex z2 = gsl_complex_rect(G-v.dat[1], v.dat[0]);
    gsl_complex z1_Y = gsl_complex_pow_real(z1, Y);
    gsl_complex z2_Y = gsl_complex_pow_real(z2, Y);

    //gsl_complex z3 = gsl_complex_rect(-0.5*sigma*sigma*v.dat[0]*v.dat[0],
    //                                  +0.5*sigma*sigma*v.dat[1]*v.dat[1]);


    gsl_complex tmp = gsl_complex_add(z1_Y,z2_Y);
    //tmp = gsl_complex_add(tmp,z3);
    tmp = gsl_complex_add_real(tmp,-std::pow(M,Y)-std::pow(G,Y));

    tmp = gsl_complex_mul_real(tmp, C*t*boost::math::tgamma(-Y));
    gsl_complex z3 = gsl_complex_rect(t*0.5*sigma*sigma*(-v.dat[0]*v.dat[0]+v.dat[1]*v.dat[1]),
                                           -t*sigma*sigma*v.dat[0]*v.dat[1]);

    tmp = gsl_complex_add(tmp, gsl_complex_rect(-v.dat[1]*t*drift,v.dat[0]*t*drift));
    tmp = gsl_complex_add(tmp, z3);
    tmp = gsl_complex_exp(tmp);

    return tmp;
}

template <typename T>
gsl_complex
CharacteristicFunction1D<T,CGMYe>::zeta_0(void)
{
    T C = processparameters.k_C;
    T G = processparameters.k_G;
    T M = processparameters.k_M;
    T Y = processparameters.k_Y;
    gsl_complex res = gsl_complex_rect(C*boost::math::tgamma(-Y)*(-Y*std::pow(M-1,Y-1)
                                       +Y*std::pow(G+1,Y-1))+drift,0.);
    return res;
}

}   // namespace lawa
