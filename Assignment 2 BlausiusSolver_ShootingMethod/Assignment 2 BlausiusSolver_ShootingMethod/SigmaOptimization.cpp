#include "SigmaOptimization.h"
#include<cmath>
#include<string>
#include<iomanip>
#include<iostream>
#include<sstream>
#include <omp.h>


SigmaOptimization::SigmaOptimization(double sigma_initial_val, double tolerance_val) : current_sigma(sigma_initial_val), tolerance(tolerance_val) {
	//Initial simulation :
	BlausiusSolver Solver(sigma_initial_val, true);
	Solver.RK4(0);
	currentfinf_sigma = Solver.GetState().f[1];
}
																			 

void SigmaOptimization::Approximation_PlusNeg() {
	DeltaSigma = coeff_deltasigma * current_sigma;
	double sigma_plus = current_sigma + 0.5 * DeltaSigma;
	double sigma_neg = current_sigma - 0.5 * DeltaSigma;

#pragma omp parallel sections
	{
	
	#pragma omp parallel section
	{
		BlausiusSolver SolverPlus(sigma_plus, false);
		SolverPlus.RK4(0); //RK4 Blausius equation solver 
		fplus_sigma = SolverPlus.GetState().f[1];
	}

	#pragma omp parallel section 
	{
		BlausiusSolver SolverNeg(sigma_neg, false);
		SolverNeg.RK4(0); //RK4 Blausius equation solver 
		fneg_sigma = SolverNeg.GetState().f[1];
	}

	}
	

	
}

void SigmaOptimization::Approximation_Sensitivity_Hessian() {
	//Approximating fplus and fneg : 
	Approximation_PlusNeg();

	//Approximating derivatives : 
	double derivative_1 = (fplus_sigma - fneg_sigma) / DeltaSigma;
	double derivative_2 = (fplus_sigma - 2 * currentfinf_sigma + fneg_sigma) / (0.25 * DeltaSigma * DeltaSigma);
	
	//Approximating sensitivity : 
	Sensitivity_L = 2 * (1 - currentfinf_sigma) * derivative_1;
	
	//Approximating Hessian : 
	Hessian_L = 2 * (derivative_2 * (1 - currentfinf_sigma) - derivative_1 * derivative_1);
}

void SigmaOptimization::UpdateValue(int max_iteration_val) {
	//Calculate initial L : 
	Objective_L = (1 - currentfinf_sigma) * (1 - currentfinf_sigma);
	int iteration = 0;
	
	//Iteration : 
	while (Objective_L > tolerance and iteration < max_iteration_val) {
		//Calculate Sensitivity and Hessian : 
		Approximation_Sensitivity_Hessian();

		//Calculate step size : 
		stepsize_sigma = -Sensitivity_L / Hessian_L;

		//Updating sigma : 
		current_sigma += stepsize_sigma;

		//Simulating again : 
		BlausiusSolver Solver(current_sigma, true);
		Solver.RK4(0);
		currentfinf_sigma = Solver.GetState().f[1];

		//Calculate objective : 
		Objective_L = (1 - currentfinf_sigma) * (1 - currentfinf_sigma);
		iteration++;
		std::cout << "At iteration " << iteration << ", sigma = " << current_sigma
			<< " and f = " << currentfinf_sigma << "\n";
	}

	//If the solution is done : 
	BlausiusSolver Solver(current_sigma, true);
	std::stringstream sigma_print;
	sigma_print << std::fixed << std::setprecision(15) << current_sigma;
	std::string filename = "Eta distribution for sigma = " + sigma_print.str() + ".csv";
	Solver.RK4(1, true, filename);
	currentfinf_sigma = Solver.GetState().f[1];
	eta99 = Solver.GetEta99();
	etastar = Solver.GetEtaStar();
	thetastar = Solver.GetThetaStar();
	ShapeFactor = Solver.GetShapeFactor();
	delta_ratio = Solver.GetDelta();
	deltastar_ratio = Solver.GetDeltaStar();
	theta_ratio = Solver.GetTheta();
	Cf_ratio = Solver.GetFrictionCoeff();
	Cd_ratio = Solver.GetDragCoeff();

	//Calculate objective : 
	Objective_L = (1 - currentfinf_sigma) * (1 - currentfinf_sigma);
	Approximation_Sensitivity_Hessian();
}

