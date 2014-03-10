
// (c) J. Poirot and P. Tankov, June-September 2006   
// Permission granted to redistribute and modify this code if the above copyright notice remains intact   
// Direct all questions concerning this code to tankov@math.jussieu.fr   
       
// Tests for the CGMY simulator   
       
//#include "cgmy.h"
#include "cgmy.hpp"
#include <fstream>   
#include <iostream>   
#include <iomanip>
//LK
#include<time.h>
#include <stdlib.h> //LK atoi
//#include <cstring>
//#include <QTime>
//#include <lawa/flensforlawa.h>
//#include <lawa/constructions/basis.h>
//#include <lawa/methods/adaptive/algorithms/localrefinement.h>
using namespace std;   

int main(int argc, char *argv[]){
	if (argc!=4) {
        cout << "Usage: " << argv[0] << " eps N1 N2" << endl;
        return 0;
    }

   
    int NrSim1   = atoi(argv[2]);
	int NrSim2   = atoi(argv[3]);	
    double eps = atof(argv[1]);

	typedef double T;
	T C1 = 1., G1 = 8, M1 = 8, Y1 = 1.25;
	T C2 = 1., G2 = 4, M2 = 4, Y2 = 1.25;

	T C3=1., G3=8, M3=8, Y3=0.2;
	T C4=1., G4=8, M4=16, Y4=1.5;

	T C5=1., G5=16, M5=8,  Y5=1.5;	
	T C6=1., G6=8., M6=16.,Y6=0.2;     
	//double eps = 0.000001; // truncation of the Levy measure   
	CGMYSimulator csim1(C1,G1,M1,Y1,eps);
	CGMYSimulator csim2(C2,G2,M2,Y2,eps);
	CGMYSimulator csim3(C3,G3,M3,Y3,eps);
	CGMYSimulator csim4(C4,G4,M4,Y4,eps);
	CGMYSimulator csim5(C5,G5,M5,Y5,eps);
    CGMYSimulator csim6(C6,G6,M6,Y6,eps);   
	//// Simulate a trajectory on the interval [0,1] with 1000 discretization points, write into a file   
	// {
	// 	double T=5;
	// 	int Nsteps = 10000;    
	// 	// double Tstep = 0.001;
	long int N= NrSim1;
	double Tstep = 1/((double) N);
	double X = 0.,X1=0.,X2=0.,X3=0.,X4=0.,X5=0.,X6=0.;   
	ofstream file("cgmySimulations.txt");
	//std::cout << "Name = " << file << std::endl;
	file << 0.0 << " "  << setw(12) << X6  << " "  << setw(12) << X1
		 << " " << setw(12) <<  X2 <<" " << setw(12) <<  X3 <<" " << setw(12) <<  X4 <<" "  << setw(12) <<  X5 << endl;   
	for(int i=1; i<N+1; i++) file << setw(1) << i*Tstep 
								  << "  " << setw(12) << (X1+=csim1.simZeroCenter(Tstep))
								  << "  " << setw(12) << (X2+=csim2.simZeroCenter(Tstep))
								  << "  " << setw(12) << (X3+=csim3.simZeroCenter(Tstep))
								  << "  " << setw(12) << (X4+=csim4.simZeroCenter(Tstep))
								  << "  " << setw(12) << (X5+=csim5.simZeroCenter(Tstep))
								  << "  " << setw(12) << (X6+=csim6.simZeroCenter(Tstep))
								 //<< "  " << setw(12) << (X6+=csim2.sim(Tstep))
								 //  << " " <<  setw(12) <<" first moment = " << csim.cumulant(1)
								  << endl;   
	// }   
	// Simulate the jump points of the process on [0,1] and write the result into a file   
	// {   
	// 	double t;   
	// 	double Tcur = 0;   
	// 	double Tmax = 1;   
	// 	double before;   
	// 	double after;   
	// 	ofstream file("cgmy_jumps.txt");   
	// 	double X=0;   
	// 	file << setw(12) << 0.0 << "  " << setw(12)<< X << endl;   
	// 	while(true){   
	// 		t = Tmax - Tcur;   
	// 		if(csim.simtojump(t,before,after)){   
	// 			X+=before;   
	// 			file << setw(12) << Tmax << "  " << setw(12)<< X << endl;   
	// 			break;   
	// 		}   
	// 		else{   
	// 			Tcur+=t;   
	// 			file << setw(12) << Tcur << "  " << setw(12)<< X+before << endl;   
	// 			file << setw(12) << Tcur << "  " << setw(12)<< X+after << endl;   
	// 			X+=after;   
	// 		}   
	// 	}   
	// } 
	
	// Compute mean and variance    
	{
		// clock_t time = clock();
			  
		// double maturity = 1;   
		// int Ntraj = NrSim2;//100000;   
		// double E=0, V=0, a;   
		// for(int i=0; i<Ntraj; i++){   
		// 	a = csim1.simZeroCenter(maturity);   
		// 	E+=(a/Ntraj);   
		// 	V+=(a*a/Ntraj);   
		// }      
		// V=(Ntraj/(Ntraj-1))*(V-E*E);   
		// //cout << "Theoretical mean: " << csim.cumulant(1)<<endl << "Empirical mean: "<< E << endl    
		// //<< "Theoretical variance: " << csim.cumulant(2)<<endl << "Empirical variance: "<< V << endl;
		// cout << "Empirical mean: " << E << endl;
		// time= clock() -time;
		// cout << "Time in seconds: " << ((float)time)/(CLOCKS_PER_SEC) << endl;


		// E=0;   
		// for(int i=0; i<Ntraj; i++){   
		// 	a = csim1.simZeroCenter(maturity);   
		// 	E+=(a/Ntraj);   
		// 	V+=(a*a/Ntraj);   
		// }      
		// V=(Ntraj/(Ntraj-1))*(V-E*E);   
		// //cout << "Theoretical mean: " << csim.cumulant(1)<<endl << "Empirical mean: "<< E << endl    
		// //<< "Theoretical variance: " << csim.cumulant(2)<<endl << "Empirical variance: "<< V << endl;
		// cout << "Sim 1 with C=" << C1 << ", G=" << G1 <<", M=" << M1 << ", Y=" << Y1 << endl
		// 	 << "Empirical mean: " << E << endl;
		// time= clock() -time;
		// cout << "Time in seconds: " << ((float)time)/(CLOCKS_PER_SEC) << endl;

		// E=0;
		
		// for(int i=0; i<Ntraj; i++){   
		// 	a = csim1.simZeroCenter(maturity);   
		// 	E+=(a/Ntraj);   
		// 	V+=(a*a/Ntraj);   
		// }      
		// V=(Ntraj/(Ntraj-1))*(V-E*E);   
		// //cout << "Theoretical mean: " << csim.cumulant(1)<<endl << "Empirical mean: "<< E << endl    
		// //<< "Theoretical variance: " << csim.cumulant(2)<<endl << "Empirical variance: "<< V << endl;
		// cout << "Sim 1 with C=" << C2 << ", G=" << G2 <<", M=" << M2 << ", Y=" << Y2 << endl
		// 	 << "Empirical mean: " << E << endl;
		// time= clock() -time;
		// cout << "Time in seconds: " << ((float)time)/(CLOCKS_PER_SEC) << endl;
	  
	}   
       
}  
