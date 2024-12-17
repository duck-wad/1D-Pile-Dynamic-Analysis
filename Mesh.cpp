#include <fstream>

#include "Mesh.h"
#include "Utils.h"

/* ELEMENT METHODS */

Element::Element(double L, double E, double A, double k, double c, double rho) {
	
	elementStiffness.resize(2, std::vector<double>(2, 0.0));
	elementMass.resize(2, std::vector<double>(2, 0.0));
	elementDamping.resize(2, std::vector<double>(2, 0.0));

	for (size_t i = 0; i < 2; i++) {
		for (size_t j = 0; j < 2; j++) {
			if (i == j) {
				elementStiffness[i][j] = (A * E / L * 1) + (k * L * 1 / 3);
				elementMass[i][j] = rho * A * L * 1 / 3;
				elementDamping[i][j] = c * L * 1 / 3;
			}
			else {
				elementStiffness[i][j] = (A * E / L * -1) + (k * L * 1 / 6);
				elementMass[i][j] = rho * A * L * 1 / 6;
				elementDamping[i][j] = c * L * 1 / 6;
			}
		}
	}

}

/* END OF ELEMENT METHOD */


/* MESH METHODS */

Mesh::Mesh(std::string fileName) {
	ReadFile(fileName);
	Discretize();
	Assemble();
}

void Mesh::ReadFile(std::string fileName) {
	std::string junk;

	std::ifstream infile(fileName);
	if (!infile) {
		std::cerr << "Error: Unable to open file." << std::endl;
	}

	infile >> junk >> length >> junk >> stiffness >> junk >> area >> junk >> spring >> junk >> damping >> junk >> density >> junk >> numelem;
	infile >> junk >> numsteps >> junk >> delta_t;

	forcedNode.resize(numsteps, 0.0);
	time.resize(numsteps, 0.0);
	for (size_t i = 0; i < numsteps; i++) {
		infile >> junk >> time[i] >> junk >> forcedNode[i];
	}
}

void Mesh::Discretize() {
	double dx = length / numelem;
	elements.reserve(numelem);
	//connectivity contains start and end position for each element
	connectivity.resize(numelem, std::vector<int>(2));
	coordinates.resize(numelem, std::vector<double>(2));

	//numelem+1 nodes
	globalStiffness.resize((numelem + 1), std::vector<double>(numelem + 1, 0.0));
	globalMass.resize((numelem + 1), std::vector<double>(numelem + 1, 0.0));
	globalDamping.resize((numelem + 1), std::vector<double>(numelem + 1, 0.0));
	globalForce.resize(numsteps, std::vector<double>(numelem+1, 0.0));
	//force is applied directly to the top node so apply directly to globalForce
	for (size_t i = 0; i < static_cast<size_t>(numsteps); i++) {
		globalForce[i][0] = forcedNode[i];
	}

	//traverse over length of pile using position
	double currentPos = 0.0;
	double startPos;
	double endPos;

	for (size_t i = 0; i < numelem; i++) {
		connectivity[i] = { static_cast<int>(i), static_cast<int>(i+1) };

		startPos = currentPos;
		endPos = currentPos + dx;
		coordinates[i] = { startPos, endPos };
		currentPos = endPos;

		elements.emplace_back(dx, stiffness, area, spring, damping, density);
	}
}

void Mesh::Assemble() {
	for (size_t i = 0; i < elements.size(); i++) {
		size_t node1 = connectivity[i][0];
		size_t node2 = connectivity[i][1];
		
		globalStiffness[node1][node1] += elements[i].getStiffness(0, 0);
		globalStiffness[node1][node2] += elements[i].getStiffness(0, 1);
		globalStiffness[node2][node1] += elements[i].getStiffness(1, 0);
		globalStiffness[node2][node2] += elements[i].getStiffness(1, 1);

		globalMass[node1][node1] += elements[i].getMass(0, 0);
		globalMass[node1][node2] += elements[i].getMass(0, 1);
		globalMass[node2][node1] += elements[i].getMass(1, 0);
		globalMass[node2][node2] += elements[i].getMass(1, 1);

		globalDamping[node1][node1] += elements[i].getDamping(0, 0);
		globalDamping[node1][node2] += elements[i].getDamping(0, 1);
		globalDamping[node2][node1] += elements[i].getDamping(1, 0);
		globalDamping[node2][node2] += elements[i].getDamping(1, 1);
	}
}


/* END OF MESH METHODS */