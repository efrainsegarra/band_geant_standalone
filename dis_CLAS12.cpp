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

double sigmainput=40.; //sigma parameter in rescattering amplitude [mb], default value
double betainput=8.; //beta parameter in rescattering amplitude [GeV^-2], default value
double epsinput=-0.5; //epsilon parameter in rescattering amplitude [], default value
double lambdainput=1.2;
double betaoffinput=8.;
//descr of off-shell rescattering:
//0=massdiff suppression (unused),
//1=dipole suppression(Jesschonek like) (unused)
//2=beta offshell parametrization (described in paper)
//3=no off-shell
//4=full off-shell
int offshellset=0;
int symm=0; //symmetric FSI in inclusive reaction

int phiavg=0; //average cross section over phi
int F_param=0; //structure functions parametrization 0=SLAC, 1=Christy&Bosted, 2=Alekhin et al leading twist

extern "C"{
    extern struct{
        char homedir[128];
        char a09file1[128];
        char a09file2[128];
    } dir_;
}//shared fortran variables


int main(int argc, char *argv[])
{
    
    //what do you want to calculate
    //0=F2*P, like Deeps data
    //1=diff cross section
    int calc=1;//atoi(argv[1]);
    int proton=1;//atoi(argv[2]); //0=DIS on neutron (spectator proton), 1=DIS on proton (spectator neutron)
    double Ein=10.9;//atof(argv[3]); //beam energy [GeV]
    double Q2=2;//atof(argv[4]); //Q^2 [GeV^2]
    double x=0.15;//atof(argv[5]); //x [] or W' [GeV] (invariant mass produced X)
    
    //some extra FSI parameters, not important for now
    int decay=0;
    int num_res=1;
    
    double pr=0.3;//atof(argv[6]); //spectator momentum [GeV]
    F_param=0;//atoi(argv[7]); //SLAC=0, Christy&bosted JLAb=1, Alekhin=2
    //non relativistic deuteron wf
    //0 - Paris wf
    //1 - AV18
    //2 - CD Bonn
    //3 - AV18*
    int which_wave=1;//atoi(argv[8]);
    offshellset=3;//atoi(argv[9]);
    betainput = 8;//atof(argv[11]);
    if(offshellset==0) epsinput=-0.5;
    if(offshellset==1) lambdainput=0;//atof(argv[12]);
    if(offshellset==2) betaoffinput=16.8;
    if(offshellset==3) epsinput=-0.5;
    if(offshellset==4) epsinput=-0.5;
    
    
    //  sanitycheck(calc,proton,F_param,which_wave,offshellset);
    
    strcpy(dir_.homedir,"/root/Desktop/PhD/EMC-SRC_exp/DeuteronDIScode"); //dir where we are running the program
    strcpy(dir_.a09file1,"/root/Desktop/PhD/EMC-SRC_exp/DeuteronDIScode");
    strcpy(dir_.a09file2,"/root/Desktop/PhD/EMC-SRC_exp/DeuteronDIScode");
    strcat(dir_.a09file1,"./grids/a09.sfs_lNNC");
    strcat(dir_.a09file2,"./grids/a09.dsfs_lNNC");
    
    double phir=180;//atof(argv[14]); //phi of spectator nucleon [degr]
    
    //input is Q^2 and x, cross section calculation [nanobarn/GeV^4]
    if(calc==1){
        
        
        // LAD angles.
        double theta_r[10] = {0};
        for(int i=0; i<=10; i++){ theta_r[i-1] = 160.5+i; }
        
        for(double Eprime = 2; Eprime<8; Eprime+=0.1){
            for(double theta_e = 5; theta_e<35; theta_e+=0.5){
                
                
                Q2 = 4*Ein*Eprime*(sin((theta_e*DEGRTORAD)/2)*sin((theta_e*DEGRTORAD)/2));
                x = Q2/(2*MASSN*(Ein-Eprime));
                double q = sqrt(Q2+(Ein-Eprime)*(Ein-Eprime));
                double theta_q = acos( (Ein-Eprime*cos(theta_e*DEGRTORAD))/q )/DEGRTORAD;
                
                for(pr=0.275+0.01/2.; pr<=0.6; pr+=0.01){
                    for(int index=0; index<10; index++){
                        double theta_pq = theta_r[index]-theta_q;
                        if (theta_pq>180)
                            theta_pq = 360 - theta_pq;
                        
                        double W_prime=sqrt( MASSD*MASSD-
                                            Q2+MASSN*MASSN+2.*MASSD*((Q2/(2*MASSN*x))-sqrt(pr*pr+MASSP*MASSP))-
                                            2.*(Q2/(2*MASSN*x))*sqrt(pr*pr+MASSP*MASSP)+
                                            2.*sqrt(Q2+(Q2/(2*MASSN*x))*(Q2/(2*MASSN*x)))*pr*cos(theta_pq*DEGRTORAD) );
                        
                        sigmainput = (25.3+53*(W_prime-MASSN))/(Q2);
                        
                        
                        if (W_prime > 2){
                            double crosstotal1 = calc_cross(Ein, Q2, x, pr, theta_pq*DEGRTORAD, phir, proton, which_wave, decay, num_res, 0);  //plane-wave
                            cout << index << " " << W_prime << " " << Q2 << " " << x << " " << pr << " " << theta_pq << " " << crosstotal1*x*Ein*(Ein-Q2/(2*MASSN*x))*pr/(PI*Q2/(2*MASSN*x)) << endl;
                        }
                        
                        
                    }
                }
            }
        }
        
        return 0; 
    }
    
}



