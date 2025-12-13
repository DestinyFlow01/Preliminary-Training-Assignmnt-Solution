#pragma once
#include <string>
#include "MainControl.h"

struct State {
    double f[3] = { 0,0,0 };
};

class BlausiusSolver {
    private : 
        State Y_val; 
        double k1[3] = { 0,0,0 };
        double k2[3] = { 0,0,0 };
        double k3[3] = { 0,0,0 };
        double k4[3] = { 0,0,0 };
        double delta_eta = 0; 
        bool Blausius_type = true;
        int num_points = 0; 

        //Boundary layer properties : 
        double eta_99 = 0;
        double BL_Params_Blausius[8] = {};
        //structure = δ*/δ, θ/δ, H, δ√Re_x / x, δ*√Re_x / x, θ√Re_x / x, C_f√Re_x, C_d√Re_x
        double sumterm_etastar = 0; 
        
    public : 
        //Constructor : 
        BlausiusSolver(double sigma, bool Blausius);

        inline void Calculate_k(double* k_calculated, int index, double* k_previous = nullptr);
        inline void RK4_OneTimeStep(int i, bool obtain99 = false);
        void RK4(bool outputCSVfile, bool obtain99 = false, std::string filename = "");

        //Accessor : 
        inline State GetState() const { return Y_val; }
        inline double GetEta99() const { return eta_99; }
        inline double GetEtaStar() const { return BL_Params_Blausius[0]; }
        inline double GetThetaStar() const { return BL_Params_Blausius[1]; }
        inline double GetShapeFactor() const { return BL_Params_Blausius[2]; }
        inline double GetDelta() const { return BL_Params_Blausius[3]; }
        inline double GetDeltaStar() const { return BL_Params_Blausius[4]; }
        inline double GetTheta() const { return BL_Params_Blausius[5]; }
        inline double GetFrictionCoeff() const { return BL_Params_Blausius[6]; }
        inline double GetDragCoeff() const { return BL_Params_Blausius[7]; }





};