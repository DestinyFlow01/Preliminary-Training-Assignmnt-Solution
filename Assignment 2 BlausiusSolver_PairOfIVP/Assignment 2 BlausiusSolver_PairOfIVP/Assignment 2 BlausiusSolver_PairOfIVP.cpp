#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <chrono>
#include<cmath>
#include "BlausiusSolver.h"
using namespace std;

double ObtainSigma();

int main()
{
    auto time_begin = chrono::high_resolution_clock::now();
    double sigma = ObtainSigma(); 

    //Solver for f :
    BlausiusSolver Solver_f(sigma, true);
    stringstream sigma_print;
    sigma_print << std::fixed << std::setprecision(15) << sigma;
    string filename = "Eta distribution for sigma = " + sigma_print.str() + ".csv";

    bool outputfilenBL = 1;
    Solver_f.RK4(outputfilenBL, true, filename);
    
    cout << fixed << setprecision(20) << "\nOptimized result : \n"
         << "sigma = " << sigma << "\n\n";
        if (outputfilenBL) {
           cout << "Boundary layer parameters : \n"
                << "eta 99% freestream = " << Solver_f.GetEta99() << "\n"
                << "etastar = " << Solver_f.GetEtaStar() << "\n"
                << "thetastar = " << Solver_f.GetThetaStar() << "\n"
                << "Shape Factor = " << Solver_f.GetShapeFactor() << "\n"
                << "delta = " << Solver_f.GetDelta() << "\n"
                << "delta star = " << Solver_f.GetDeltaStar() << "\n"
                << "theta = " << Solver_f.GetTheta() << "\n"
                << "Friction coefficient = " << Solver_f.GetFrictionCoeff() << "\n"
                << "Drag coefficient = " << Solver_f.GetDragCoeff() << "\n";
        }
    
    auto time_end = chrono::high_resolution_clock::now();
    double time_elapsed_ms = chrono::duration_cast<chrono::milliseconds>(time_end - time_begin).count();

    //At the end : 
    double sigma_reference = 0.469599988361013304509;
    double relative_error = abs(sigma_reference - sigma) / sigma_reference;
    cout << "\nThe relative error compared to reference value " << sigma_reference << setprecision(10) << " is " << relative_error * 100 << "%\n";


    cout << "\n\nTime needed for calculation = " << fixed << setprecision(3) << time_elapsed_ms << " ms\n";

    char ch; cin >> ch;
}

double ObtainSigma()
{
    //Solver for F :
    BlausiusSolver Solver_F(1, false);
    Solver_F.RK4(0);
    double sigma = pow(Solver_F.GetState().f[1], -1.5);
    return sigma;
}
