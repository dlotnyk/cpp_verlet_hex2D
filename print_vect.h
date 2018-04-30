#ifndef PRINT_VECT_H_
#define PRINT_VECT_H_ 
#include <vector>
#include <iostream>
#include <string>
#include <stdio.h>

struct point
{
	std::vector<double> x; // cartesian
	std::vector<double> y;
	std::vector<int> a; // hexagonal
	std::vector<int> b;
	void initialization(int num) {
		a = std::vector<int>(num, 0);
		b = std::vector<int>(num, 0);
		x = std::vector<double>(num, 0);
		y = std::vector<double>(num, 0);
	}
	void clear() {
		x.clear();
		y.clear();
		a.clear();
		b.clear();
	}
	void delete_n(int n) {
		x.erase(x.begin() + n);
		y.erase(y.begin() + n);
		a.erase(a.begin() + n);
		b.erase(b.begin() + n);
}
};

class print_vect
{
public:
	void print_itself(const std::vector<double> &);
	void print_itself(const std::vector<double> &, const std::vector<double> &, std::string);
	void print_itself(const std::vector<int> &);
	void print_itself(point);
	void print_itself(point, std::string);
};

#endif PRINT_VECT_H_