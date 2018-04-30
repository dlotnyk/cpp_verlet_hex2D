#ifndef LATTICE_GRID_H_
#define LATTICE_GRID_H_

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <new>
#include <iomanip>
#include <chrono>
#include "print_vect.h"

struct dot {
	double x, y;
};

class lattice_grid: public print_vect
{
private:
	point base_hex_lattice; // 1-st coordinate circle in hex coordinates
	int number_p;
protected:
	dot cent;
	point neibor, lattice, sample;
public:
	
	void hex_transform(point &);
	void gen_lattice(); // generate lattice
	void gen_sample();
	void neighbor(int ini_point);

	lattice_grid();
	lattice_grid(const int &nn);
	~lattice_grid();
};

#endif LATTICE_GRID_H_