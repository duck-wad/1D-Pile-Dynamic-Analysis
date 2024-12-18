#include <iostream>

#include "Mesh.h"

int main() {

	Mesh mesh("INPUT.txt");
	mesh.Solve();
	mesh.PrintResults();

	return 0;
}