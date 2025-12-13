#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <chrono>
#include<omp.h>
#include "BlausiusSolver.h"
#include "SigmaOptimization.h"
#include "MainControl.h"
using namespace std;

int main()
{  
    omp_set_num_threads(num_threads);
    SigmaOptimization Optimization(sigma_initial, tolerance); //This is Blausius connection parameter, which is the second derivative of f
    auto time_begin = chrono::high_resolution_clock::now();
    Optimization.UpdateValue(max_iteration);
    auto time_end = chrono::high_resolution_clock::now();
    double time_elapsed_ms = chrono::duration_cast<chrono::milliseconds>(time_end - time_begin).count();

    //At the end : 
    cout << "Time needed for calculation = " << fixed << setprecision(3) << time_elapsed_ms << " ms\n";
    
    cout << fixed << setprecision(20) << "\nOptimized result : \n"
        << "sigma = " << Optimization.Get_currentsigma() << "\n"
        << "Current f = " << Optimization.Get_f() << "\n"
        << "Objective function = " << Optimization.Get_Objective() << "\n"
        << "Sensitivity = " << Optimization.Get_Sensitivity() << "\n"
        << "Hessian = " << Optimization.Get_Hessian() << "\n\n"
        << "Boundary layer parameters : \n"
        << "eta 99% freestream = " << Optimization.GetEta99() << "\n"
        << "etastar = " << Optimization.GetEtaStar() << "\n"
        << "thetastar = " << Optimization.GetThetaStar() << "\n"
        << "Shape Factor = " << Optimization.GetShapeFactor() << "\n"
        << "delta = " << Optimization.GetDelta() << "\n"
        << "delta star = " << Optimization.GetDeltaStar() << "\n"
        << "theta = " << Optimization.GetTheta() << "\n"
        << "Friction coefficient = " << Optimization.GetFrictionCoeff() << "\n"
        << "Drag coefficient = " << Optimization.GetDragCoeff() << "\n";
        


    double sigma_reference = 0.469599988361013304509; 
    double relative_error = abs(sigma_reference - Optimization.Get_currentsigma()) / sigma_reference; 
    cout << "\nThe relative error compared to reference value " << sigma_reference <<setprecision(10) << " is " << relative_error * 100 << "%\n";

    //char ch; cin >> ch;
}
