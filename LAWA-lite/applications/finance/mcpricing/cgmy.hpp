// (c) J. Poirot and P. Tankov, June-September 2006 
// Permission granted to redistribute and modify this code if the above copyright notice remains intact 
// Direct all questions concerning this code to tankov@math.jussieu.fr 
 
 
// Simulation of the CGMY process with Levy measure truncated at level eps  
// (jumps smaller than eps in absolute value are replaced with their mean) 
// Uses the algorithm in Madan and Yor (), see also Poirot and Tankov (2006) 
// The Levy measure of CGMY process is  
// C exp(-M x) / x^{1+Y} 1_{x>0} + C exp(-G |x|) / |x|^{1+Y} 1_{x<0} 
// the gamma parameter of the Levy triplet is chosen in such way that 
// the mean of X_1 is equal to C gamma(1-Y) (M^{Y-1}-G^{Y-1}) 
// this is the natural 'zero drift' version of the process  
// corresponding to using a subordinator without drift 



#ifndef _CGMYSIM 
#define _CGMYSIM

//LK from cgmy.tcc
#include <boost/math/special_functions/gamma.hpp>
#include <cmath>   
#include <cstdlib> //LK, atoi, srand, rand
#include <time.h> //LK, time
//LK, different gamma function



//-----
 
class CGMYSimulator{ 
	const double C, G, M, Y; 
	const double eps; 
	double A, B, d, lambda, P; 
public: 
	CGMYSimulator(double _C, double _G, double _M, double _Y, double _eps); 
	long double sim(double t); // simulate a t-increment of the truncated CGMY process
	long double simZeroCenter(double t);  //LK: simulates a t-increment of the truncated CGMY process with Zero center!
	
	// returns the increment value 
	bool simtojump(double & t, double & before, double & after); // simulate the (truncated) CGMY process 
	// up to the first jump or up to time t, if the jump arrives after t 
	// on entry, t contains the time step 
	// returns true if the jump arrives after t and false otherwise 
	// if true is returned, before contains the increment value 
	// if false is returned, t contains the jump moment, before contains the process value before the jump 
	// and after contains the process value after the jump 
	long double cumulant(int n); 
	// returns the n-th cumulant of X_1 (computed theoretically): 
	// cumulant(1) corresponds to the mean, cumulant(2) to the variance etc.
	

	//LK: implementations of cgmy.cpp
private:
	static const double Pi;
	long double theoreticalMean;
	double Uniform ()   
	{
		return (double) rand() /RAND_MAX;   
	}   
       
	double StandardGaussian()   
	{   
		static bool haveVar = false;   
		static double var;   
		if(haveVar) {   
			haveVar = false;   
			return var;   
		}   
		double u, v;   
		u=Uniform();   
		v=Uniform();   
		var=sqrt(-2*log(u))*cos(2*Pi*v);   
		haveVar = true;   
		return sqrt(-2*log(u))*sin(2*Pi*v);   
	}   
       
       
	double gammLN(double xx) {   
		// Returns lnG(xx)   
		// Lanczos formula from Numerical Recipes in C   
		double x,y,tmp,ser;   
		const double coef[6]={76.18009172947146,-86.50532032941677,   
							  24.01409824083091,-1.231739572450155,   
							  0.1208650973866179e-2,-0.5395239384953e-5};   
		int j;   
		y=x=xx;   
		tmp=x+5.5;   
		tmp -= (x+0.5)*std::log(tmp);   
		ser=1.000000000190015;   
		for (j=0;j<=5;j++) ser += coef[j]/++y;   
		return -tmp+std::log(2.5066282746310005*ser/x);   
	}   
       
	inline double implemented_gamma(double x) {return std::exp(gammLN(x));} // gamma function   
       
	class InfExcept{};   
       
	double DegHypergeometric1(double a, double b, double z)   
	{   
		int i=0;   
		double H=1;   
		double c=1;   
		while(fabs(c)>0.00000001)   
			{   
				c *= ((a+i)/(b+i))*z/(i+1);   
				if(fabs(c)>1e20) throw InfExcept();   
				H+=c;   
				i++;   
				if(i>300) {   
					//          std::cout << "Maximum iteration number reached in DegHypergeometric with z="<<z<<"; H="<<H<<std::endl;   
					break;   
				}   
			}   
       
		return H;   
	}   
       
	double ParCylFunction1 (double p, double z)   
	{     
             
		int i;   
		double u,v;   
		double S=1, c=1;   
		double g=implemented_gamma(-p/2+0.5);   
		double h=implemented_gamma(-p/2);   
		double d;   
		if(z<40)   
			{   
				u= DegHypergeometric1(-p/2,0.5,z*z/2);   
				v= DegHypergeometric1(0.5-p/2,1.5,z*z/2);   
				d=pow(2,p/2)*exp(-z*z/4)*((sqrt(Pi)*u/g)-sqrt(2*Pi)*z*v/h);   
				return d;   
			}   
		else{   
			for (i=1;i<20;i++)   
				{    
					c*= -(p-i+1)*(p-i)/(2*i*pow(z,2*i));   
					S+=c;}   
                
			return exp(-z*z/4)*pow(z,p)*S;}   
            
	}   
       
       
            
	double IntegralI(double Y,double a, double lambda)   
	{   
		double val;   
		val= pow(2*lambda, (-Y)/2)*implemented_gamma(Y)*exp(a*a/(8*lambda))*ParCylFunction1(-Y, (a/sqrt(2*lambda)));   
           
		return val;   
	}   
       
	double h(double y, double Y, double A,double B)   
	{   
		double val1;   
		val1=(exp((A*A-B*B)*y/2)*implemented_gamma((Y+1)/2)*pow(2,Y)*pow(B*B*y/2,Y/2)*IntegralI(Y,B*B*y,B*B*y/2))/(implemented_gamma(Y)*implemented_gamma(0.5));   
		return val1;   
	}


};







           ///////////////////////////////////////////////////
           //          Class Implementations                //
           ///////////////////////////////////////////////////


const double CGMYSimulator::Pi=3.141592653589793;   

    CGMYSimulator::CGMYSimulator(double _C, double _G, double _M, double _Y, double _eps) :    
    C(_C), G(_G), M(_M), Y(_Y), eps(_eps){
		
        P = implemented_gamma(0.5)/(pow(2.0,Y/2)*implemented_gamma((Y+1)/2));    
        A = (G-M)/2;   
        B = (G+M)/2;   
        d = P*C*pow(eps,1-Y/2)/(1-Y/2);   
        lambda = 2*C*P/(Y*pow(eps,Y/2));
		srand(time(NULL));
		theoreticalMean = this->cumulant(1);
    }   
       
    long double CGMYSimulator::sim(double t){   
        double u1, u2, u3, y;   
        double H=d*t;   
        double K=0;   
        while(K<t)   
            {      
                u1=Uniform();   
                u2=Uniform();   
                u3=Uniform();   
                K-=log (u2)/lambda;   
                if(K>t) break;   
                y= eps/(pow(u1,2/Y));   
                try{   
                    if (h(y,Y,A,B)>u3) H+=y;   
                }   
                catch(InfExcept){   
                }   
            }   
      return A*H+sqrt(H)*StandardGaussian();       
    }

    //LK: zero center simulation!
    long double CGMYSimulator::simZeroCenter(double t){
		return this->sim(t) - theoreticalMean*t;
    }
       
    bool CGMYSimulator::simtojump(double & t, double & before, double & after){   
        double y;   
        double K=0;   
        while(K<t)   
                {      
                     double u1=Uniform();   
                     double u2=Uniform();   
                     double u3=Uniform();   
                     double f=-log (u2)/lambda;   
                     K += f;   
                     if(K>t) break;   
                     y= eps/(pow(u1,2/Y));   
                     try{   
                         if (h(y,Y,A,B)>u3){   
                            t = K;   
                            double sg1 = StandardGaussian();   
                            double sg2 = StandardGaussian();   
                            before = A*d*t+sqrt(d*t)*sg1;   
                            after = before + A*y+sqrt(y)*sg2;   
                            return false;   
                         }   
                     }   
                     catch(InfExcept){   
                     }   
                           
                }   
        double sg = StandardGaussian();   
        after = before = A*d*t+sqrt(d*t)*sg;   
        return true;   
    }   
       
    long double CGMYSimulator::cumulant(int n){
		// LK, test, was: implemented_gamma, is boost::math::tgamma
        if(n%2 == 0){
			return C*boost::math::tgamma(-Y+n)*(pow(M,Y-n) + pow(G,Y-n));
		}
		else{

			// std::cout << " tgamma(1-Y) = " << boost::math::tgamma(-Y+n)
			// 		  << " M^{Y-1} = " << pow(M,Y-n)
			// 		  << " G^{Y-1} = " << pow(G,Y-n)
			// 		  << " E = " <<boost::math::tgamma(-Y+n)*(pow(M,Y-n) - pow(G,Y-n))   //C*implemented_gamma(-Y+n)*(pow(M,Y-n) - pow(G,Y-n)) 
			//		  << " Theoretical Mean = " << theoreticalMean << std::endl;
			//std::cout << " cumulant: " << C*implemented_gamma(-Y+n)*(pow(M,Y-n) - pow(G,Y-n)) << std::endl ;
			return C*boost::math::tgamma(-Y+n)*(pow(M,Y-n) - pow(G,Y-n));
		}   
	}   




#endif 
