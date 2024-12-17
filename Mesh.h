#pragma once

#include <iostream>
#include <vector>

class Element {
public:
	Element(double L, double E, double A, double k, double c, double rho);
	//return row i column j of matrix
	double getStiffness(size_t i, size_t j) { return elementStiffness[i][j]; }
	double getMass(size_t i, size_t j) { return elementMass[i][j]; }
	double getDamping(size_t i, size_t j) { return elementDamping[i][j]; }

protected:
	std::vector<std::vector<double>> elementStiffness;
	std::vector<std::vector<double>> elementMass;
	std::vector<std::vector<double>> elementDamping;
};

class Mesh {
public: 
	Mesh(std::string fileName);

	void ReadFile(std::string fileName);
	void Discretize();
	void Assemble();
	void ApplyBC();
	void Solve();
	void PrintResults();

protected:
	//connectivity matrix is 1:1 since 1D 
	std::vector<std::vector<int>> connectivity;
	std::vector<std::vector<double>> coordinates;
	int numelem;
	//pile parameters
	double length;
	double stiffness;
	double area;
	double spring;
	double damping;
	double density;
	//time parameters
	int numsteps;
	double delta_t;
	std::vector<double> time;
	std::vector<double> forcedNode;
	//each element object will contain its own stiffness, mass, damping matrices, and force vector
	std::vector<Element> elements;
	//global stiffness, mass, damping matrices
	std::vector<std::vector<double>> globalStiffness;
	std::vector<std::vector<double>> globalMass;
	std::vector<std::vector<double>> globalDamping;
	//global force vector is a vector of vectors, since force is time dependent
	std::vector<std::vector<double>> globalForce;
	//displacement, velocity, and acceleration will have same dimensions as globalForce
	std::vector<std::vector<double>> displacement;
	std::vector<std::vector<double>> velocity;
	std::vector<std::vector<double>> acceleration;
};