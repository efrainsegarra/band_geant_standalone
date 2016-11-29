#include <iostream>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <complex>
#include <fstream>
#include <cstdarg>
#include <time.h>

using namespace std;

#include "constants.h"
#include "deuteronwf.h"
#include "crossdis.h"
#include "crossincl.h"
double sigmainput=40.;
double betainput=8.;
double epsinput=-0.5;
double lambdainput=1.2;
double betaoffinput=8.;
int offshellset=0;

int phiavg=1;
int F_param=2;
int symm=0;
extern "C"{
    extern struct{
        char homedir[128];
        char a09file1[128];
        char a09file2[128];
    } dir_;
}//shared fortran variables

int main(int argc, char *argv[])
{
    
    
    int calc=0;//atoi(argv[1]);
    int which_wave=1;//atoi(argv[6]);
    F_param=0;//atoi(argv[5]);
    int decay=0;
    int num_res=1;
    int offshellset=3;
    sanitycheck(calc,0,F_param,which_wave,offshellset);
    
    double *c, *d, *m;
    int c_length=0;  //parametriz array dim
    get_wf_param(&c, &d, &m,c_length, which_wave);
    
    double Ein=10.9;//atof(argv[2]); //incoming beam energy [GeV]
    double Q2=2;//atof(argv[3]); //[GeV^2]
    double x=0.1;//atof(argv[4]);
    strcpy(dir_.homedir,"/root/Desktop/PhD/EMC-SRC_exp/DeuteronDIScode"/*argv[7]*/); //dir where we are running the program
    strcpy(dir_.a09file1,"/root/Desktop/PhD/EMC-SRC_exp/DeuteronDIScode"/*argv[7]*/);
    strcpy(dir_.a09file2,"/root/Desktop/PhD/EMC-SRC_exp/DeuteronDIScode"/*argv[7]*/);
    strcat(dir_.a09file1,"/grids/a09.sfs_lNNC");
    strcat(dir_.a09file2,"/grids/a09.dsfs_lNNC");
    
    
    double QEpw, DISpw, DISfsi, DISfsi2;
    
    
    
    
    
    double Eprime_central = 6.8;
    double Eprime_bin = 0.1;
    
    double theta_e_central = 14;
    double theta_e_bin = 0.2;
    double theta_e_acceptance = 1.6;
    
    for(double Eprime = 2; Eprime<8; Eprime+=0.1){
        
        for(double theta_e = 5; theta_e<35; theta_e+=0.5){
            
            //      double Eprime = 4.5; double theta_e = 17;
            Q2 = 4*Ein*Eprime*(sin((theta_e*DEGRTORAD)/2)*sin((theta_e*DEGRTORAD)/2));
            x = Q2/(2*MASSN*(Ein-Eprime));
            
            //calc inclusive
            calc_inclusive2(QEpw, DISpw, DISfsi, DISfsi2, Ein, Q2,x, c, d, m, c_length, which_wave, offshellset, decay, num_res, calc);
            double dsigma=DISpw+QEpw-DISfsi;
            cout << Q2 << " " << x << " " /*<< QEpw*x*Ein*(Ein-Q2/(2*MASSN*x))/(PI*Q2/(2*MASSN*x)) << " " << DISpw*x*Ein*(Ein-Q2/(2*MASSN*x))/(PI*Q2/(2*MASSN*x)) << " "*/  << dsigma*x*Ein*(Ein-Q2/(2*MASSN*x))/(PI*Q2/(2*MASSN*x)) << endl;
            
            
        }
    }

    
    
    
    
    delete [] c; delete [] d; delete [] m;
    
    
}
