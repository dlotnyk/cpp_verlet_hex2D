#ifndef COORD_H_
#define COORD_H_

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <new>
#include <iomanip>
#include <chrono>
//#include "print_vect.h"
#include "lattice_grid.h"
#include "timer.h"
struct ini_condit 
{
	double k;
	double m;
	double dt;
	double x0;
	double y0;
	double eta;
};

struct values
{
	std::vector<double> x;
	std::vector<double> y;
	std::vector<double> vx;
	std::vector<double> vy;
	std::vector<double> ax;
	std::vector<double> ay;
	void initialization(int num) {
		x = std::vector<double>(num, 0);
		y = std::vector<double>(num, 0);
		vx = std::vector<double>(num, 0);
		vy = std::vector<double>(num, 0);
		ax = std::vector<double>(num, 0);
		ay = std::vector<double>(num, 0);
	}
	void clear() {
		x.clear();
		y.clear();
		vx.clear();
		vy.clear();
		ax.clear();
		ay.clear();
	}
};

class coord: public lattice_grid, public timer
{
private:
	//lattice_grid lat;
	int number_p;
	point acceler, velo;
protected:
	point pp, vor_coor;
	std::vector<double> ampl;
	std::vector<double> orient;
	ini_condit beg;

public:
	values val;
	coord();
	coord(const int &nn);
	coord(const int &nn, const double &m, const double &k, const double &dt, const double &x0, const double &y0, const double &eta);
	coord(const int &nn, const double &m, const double &k, const double &dt, std::vector<double> &x0, std::vector<double> &y0, const double &eta);
	~coord();
	//void print_itself(point);
	void set_values(const std::vector<double> &, const std::vector<double> &);
	std::vector<double> amplitude_point(const int &);
	void change_lattice(double x, double y); // change only center
	void change_lattice(double x, double y, int a); // change all in one direction
	void change_lattice(); //make random displacement
	void change_lattice(std::vector<double> &x, std::vector<double> &y); // if we want to asign initial lattice
	void change_velo(); // default change in algorithm
	void velo_begin(double vx, double vy); // change velosity of center vortex
	void velo_begin(double vx, double vy, int a); // change velosity of all vorticies
	void change_coord(); // change coord in algorithm
	void accelaration(); // calculation of accelatation according to displacements
	void verlet(); // verlet algorithm
	void get_values(int num); // get coordinates, velocita, accel for num - number of vorticies
	
};

#endif COORD_H_