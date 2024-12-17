#pragma once

#include <iostream>
#include <vector>

//use average acceleration method for time integration 
void AverageAccelerationMethod(std::vector<std::vector<double>>& D, std::vector<std::vector<double>>& V, std::vector<std::vector<double>>& A, const std::vector<std::vector<double>>& F, const std::vector<double>& D_i, const std::vector<double>& V_i, const int numsteps, const double delta_t, const std::vector<std::vector<double>>& M, const std::vector<std::vector<double>>& C, const std::vector<std::vector<double>>& K, const int nodes);
