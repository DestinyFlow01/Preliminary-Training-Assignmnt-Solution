#pragma once
#include "BlausiusSolver.h"
#include "MainControl.h"

class SigmaOptimization {
	private : 
		//Current value : 
		double currentfinf_sigma = 0;
		double current_sigma = sigma_initial;
		double stepsize_sigma = 0.1;
		double tolerance = 1e-5;
		double DeltaSigma = 1; 

		//Approximation purpose : 
		double fplus_sigma = 0; 
		double fneg_sigma = 0;

		//Objective function : 
		double Objective_L = 1; 
		double Sensitivity_L = 1; 
		double Hessian_L = 1;

		//Solver parameter : 
		double eta99 = 0; 
		//structure = δ*/δ, θ/δ, H, δ√Re_x / x, δ*√Re_x / x, θ√Re_x / x, C_f√Re_x, C_d√Re_x
		double etastar = 0; 
		double thetastar = 0; 
		double ShapeFactor = 0; 
		double delta_ratio  = 0; 
		double deltastar_ratio = 0; 
		double theta_ratio = 0; 
		double Cf_ratio = 0; 
		double Cd_ratio = 0; 
	public : 
		SigmaOptimization(double sigma_initial_val, double tolerance_val);

		//Getter to external (optimization parameters): 
		double Get_f() const { return currentfinf_sigma; };
		double Get_Objective() const { return Objective_L; };
		double Get_Sensitivity() const { return Sensitivity_L; };
		double Get_Hessian() const { return Hessian_L; };
		double Get_currentsigma() const { return current_sigma; };

		//Getter for external (BL Parameters)
		inline double GetEta99() const { return eta99; }
		inline double GetEtaStar() const { return etastar; }
		inline double GetThetaStar() const { return thetastar; }
		inline double GetShapeFactor() const { return ShapeFactor; }
		inline double GetDelta() const { return delta_ratio; }
		inline double GetDeltaStar() const { return deltastar_ratio; }
		inline double GetTheta() const { return theta_ratio; }
		inline double GetFrictionCoeff() const { return Cf_ratio; }
		inline double GetDragCoeff() const { return Cd_ratio; }

		//Approximation function 
		void Approximation_PlusNeg();
		void Approximation_Sensitivity_Hessian();
		void UpdateValue(int max_iteration_val);
};

