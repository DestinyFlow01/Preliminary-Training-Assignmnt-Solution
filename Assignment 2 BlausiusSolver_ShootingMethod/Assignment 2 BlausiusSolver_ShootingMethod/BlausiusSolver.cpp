#include <iostream>
#include <iomanip>
#include<string>
#include <sstream>
#include<fstream>
#include "BlausiusSolver.h"
#include "SigmaOptimization.h"

using namespace std;

//Member function definition for class BlausiusSolver

//Parametric Constructor 
BlausiusSolver::BlausiusSolver(double sigma, bool Blausius) : Y_val( {0,0,sigma} ), Blausius_type(Blausius) {
    delta_eta = (Blausius) ? delta_eta_Blausius : delta_eta_ApproxDeriv;
    num_points = (Blausius) ? num_points_Blausius : num_points_ApproxDeriv;

    if (Blausius) {
        BL_Params_Blausius[1] = sigma;
        BL_Params_Blausius[5] = sigma * sqrt(2);
        BL_Params_Blausius[6] = BL_Params_Blausius[5];
        BL_Params_Blausius[7] = 2 * BL_Params_Blausius[6];
    }
        
}

void BlausiusSolver::Calculate_k(double* k_calculated, int index, double* k_previous) {
    /*
    This is practically calculating F, given the input value and previous k value.
    */
    switch (index) {
    case 1: {
        k_calculated[0] = Y_val.f[1];
        k_calculated[1] = Y_val.f[2];
        k_calculated[2] = -Y_val.f[0] * Y_val.f[2];
        break;
    }

    case 2:
    case 3: {
        k_calculated[0] = Y_val.f[1] + 0.5 * delta_eta * k_previous[1];
        k_calculated[1] = Y_val.f[2] + 0.5 * delta_eta * k_previous[2];
        k_calculated[2] = -(Y_val.f[0] + 0.5 * delta_eta * k_previous[0]) * (Y_val.f[2] + 0.5 * delta_eta * k_previous[2]);
        break;
    }

    case 4: {
        k_calculated[0] = Y_val.f[1] + delta_eta * k_previous[1];
        k_calculated[1] = Y_val.f[2] + delta_eta * k_previous[2];
        k_calculated[2] = -(Y_val.f[0] + delta_eta * k_previous[0]) * (Y_val.f[2] + delta_eta * k_previous[2]);
        break;
    }
    default: {
        cerr << "Index unidentified"; break;
    }

    }

}

inline void BlausiusSolver::RK4_OneTimeStep(int i, bool obtain99nBL) {
    Calculate_k(k1, 1);      //Calculate k1
    Calculate_k(k2, 2, k1);  //Calculate k2
    Calculate_k(k3, 3, k2);  //Calculate k3
    Calculate_k(k4, 4, k3);  //Calculate k4

    double Yval_f1_previous = Y_val.f[1];
    bool check = true;
    //Runge kutta updating mechanism 
    #pragma omp parallel for 
    for (int idx = 0; idx < 3; idx++) Y_val.f[idx] += delta_eta * (k1[idx] + 2 * k2[idx] + 2 * k3[idx] + k4[idx]) / 6.;

    if (obtain99nBL) {
        if (check) {
            if (Y_val.f[1] > 0.99 and Yval_f1_previous < 0.99) {
                eta_99 = (i - 1) * delta_eta + (0.99 - Yval_f1_previous) * delta_eta / (Y_val.f[1] - Yval_f1_previous);
                check = false;
            }
        }

        sumterm_etastar += 0.5 * ((1 - Y_val.f[1]) + (1 - Yval_f1_previous));
    }
    
}   

//Solving Blausius equation and outputting CSV
void BlausiusSolver::RK4(bool outputCSVfile, bool obtain99nBL, string filename) {
    if (outputCSVfile) {
        ofstream ofs;
        ofs.open(filename, ios::out);
        ofs << "eta,f,f',f''\n";
        ofs << 0 << "," << Y_val.f[0] << "," << Y_val.f[1] << "," << Y_val.f[2] << "\n";
        ofs << fixed << setprecision(6);

        for (int i = 1; i < num_points; i++) {
            RK4_OneTimeStep(i, obtain99nBL);
            ofs << i * delta_eta << "," << Y_val.f[0] << "," << Y_val.f[1] << "," << Y_val.f[2] << "\n";
            if (abs(1 - Y_val.f[1]) < 1e-12) break;
        }
        ofs.close();
    }
    else {
        for (int i = 1; i < num_points; i++)
        {
            RK4_OneTimeStep(i, obtain99nBL);
            //if (abs(1 - Y_val.f[1]) < 1e-12) break;
        }
    }   

    if (obtain99nBL) {
        BL_Params_Blausius[0] = sumterm_etastar * delta_eta;
        BL_Params_Blausius[2] = BL_Params_Blausius[0] / BL_Params_Blausius[1];
        BL_Params_Blausius[3] = eta_99 * sqrt(2);
        BL_Params_Blausius[4] = sqrt(2) * BL_Params_Blausius[0];
    }
    
}