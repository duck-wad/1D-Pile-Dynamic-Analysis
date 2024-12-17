#include <cmath>

#include "NewmarkBeta.h"
#include "Utils.h"

void AverageAccelerationMethod(std::vector<std::vector<double>>& D, std::vector<std::vector<double>>& V, std::vector<std::vector<double>>& A, const std::vector<std::vector<double>>& F, const std::vector<double>& D_i, const std::vector<double>& V_i, const int numsteps, const double delta_t, const std::vector<std::vector<double>>& M, const std::vector<std::vector<double>>& C, const std::vector<std::vector<double>>& K, const int nodes) {

	if (F.size() != numsteps) {
		throw std::invalid_argument("Forcing function needs to be discretized to same number of points as numsteps");
	}

	D.assign(numsteps, std::vector<double>(nodes, 0.0));
	V.assign(numsteps, std::vector<double>(nodes, 0.0));
	A.assign(numsteps, std::vector<double>(nodes, 0.0));

	std::vector<std::vector<double>> K_ = M * (4.0 / (delta_t * delta_t)) + C * (2.0 / delta_t) + K;
	std::vector<double> P_;

	for (size_t i = 0; i < static_cast<size_t>(numsteps); i++) {
		if (i == 0) {
			D[i] = D_i;
			V[i] = V_i;
			A[i] = invertMatrix(M) * (F[i] - K * D[i] - C * V[i]);
		}
		else {
			P_ = F[i] + M * (D[i - 1] * (4.0 / (delta_t * delta_t)) + V[i - 1] * (4.0 / delta_t) + A[i - 1])
				+ C * (D[i - 1] * (2.0 / delta_t) + V[i - 1]);
			D[i] = invertMatrix(K_) * P_;
			V[i] = (D[i] - D[i - 1]) * (2.0 / delta_t) - V[i - 1];
			A[i] = (D[i] - D[i - 1] - V[i - 1] * delta_t) * (4.0 / (delta_t * delta_t)) - A[i - 1];
		}
	}
}