// (c) J. Poirot and P. Tankov, June-September 2006   
// Permission granted to redistribute and modify this code if the above copyright notice remains intact   
// Direct all questions concerning this code to tankov@math.jussieu.fr   
       
// Tests for the CGMY simulator   
       
#include "cgmy.h"   
#include <fstream>   
#include <iostream>   
#include <iomanip>   
using namespace std;   

int main(){   
	double C=0.5;   
	double G=3;   
	double M=5;   
	double Y=0.5;     
	double eps = 0.000001; // truncation of the Levy measure   
	CGMYSimulator csim(C,G,M,Y,eps);   
       
	// Simulate a trajectory on the interval [0,1] with 1000 discretization points, write into a file   
	{   
		int Nsteps = 1000;    
		double Tstep = 0.001;   
		double X = 0;   
		ofstream file("cgmy.txt");   
		file << setw(12) << 0.0 << "  " << setw(12)<< X << endl;   
		for(int i=1; i<Nsteps+1; i++) file << setw(12) << i*Tstep << "  " << setw(12)<< (X+=csim.sim(Tstep)) << endl;   
	}   
	// Simulate the jump points of the process on [0,1] and write the result into a file   
	{   
		double t;   
		double Tcur = 0;   
		double Tmax = 1;   
		double before;   
		double after;   
		ofstream file("cgmy_jumps.txt");   
		double X=0;   
		file << setw(12) << 0.0 << "  " << setw(12)<< X << endl;   
		while(true){   
			t = Tmax - Tcur;   
			if(csim.simtojump(t,before,after)){   
				X+=before;   
				file << setw(12) << Tmax << "  " << setw(12)<< X << endl;   
				break;   
			}   
			else{   
				Tcur+=t;   
				file << setw(12) << Tcur << "  " << setw(12)<< X+before << endl;   
				file << setw(12) << Tcur << "  " << setw(12)<< X+after << endl;   
				X+=after;   
			}   
		}   
	}   
	// Compute mean and variance    
	{   
		double T = 1;   
		int Ntraj = 10000;   
		double E=0, V=0, a;   
		for(int i=0; i<Ntraj; i++){   
			a = csim.sim(T);   
			E+=(a/Ntraj);   
			V+=(a*a/Ntraj);   
		}      
		V=(Ntraj/(Ntraj-1))*(V-E*E);   
		cout << "Theoretical mean: " << csim.cumulant(1)<<endl << "Empirical mean: "<< E << endl    
			 << "Theoretical variance: " << csim.cumulant(2)<<endl << "Empirical variance: "<< V << endl;   
	}   
       
}  
