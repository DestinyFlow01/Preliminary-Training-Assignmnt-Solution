#pragma once

//Constants for Blausius : 
static constexpr int num_points_Blausius = 1e4 + 1; //This one is for the main simulation, not used for approximation
static constexpr double eta_max = 50;
static constexpr double delta_eta_Blausius = eta_max / (num_points_Blausius - 1);

//Constants for optimization : 
static constexpr double sigma_initial = 0.3;
static constexpr int num_points_ApproxDeriv = 1e3 + 1;
static constexpr double delta_eta_ApproxDeriv = eta_max / (num_points_ApproxDeriv - 1);
static constexpr double tolerance = 1e-250;
static constexpr int max_iteration = 1e5;
static constexpr double coeff_deltasigma = 0.001; 

//Parallelization : 
static constexpr int num_threads = 2; 