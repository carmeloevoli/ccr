#include <iostream>

// #include "TEvolve.h"
#include "utilities.h"

using namespace std;

int main() {

	print_timescales("output/timescales_at_z10.txt", 10);

	print_timescales("output/timescales_at_z20.txt", 20);

/*	TEvolve* E = new TEvolve(400,"var");

	E->init_energy_grid();
	E->init_equation_vectors();
	E->init_time_scales_vectors();
	E->print_time_scales();

	E->evolve();

	delete E;
*/
	return 0;
}
