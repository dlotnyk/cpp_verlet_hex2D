#include "print_vect.h"


void print_vect::print_itself(const std::vector<double> &x) {
		for (double n : x)
			std::cout << "X = " << n << '\n';
		printf("Size is %d \n", x.size());
}
void print_vect::print_itself(const std::vector<int> &x) {
	for (int n : x)
		std::cout << "A = " << n << '\n';
	printf("Size is %d \n", x.size());
}

void print_vect::print_itself(point pa)
{
	for (int ii = 0; ii < pa.a.size(); ii++)
	{
		printf("a = %4d \t b = %4d \t x = %4.4f \t y = %4.4f \n", pa.a[ii], pa.b[ii], pa.x[ii], pa.y[ii]);
	}
	printf("size is %d \n", pa.x.size());
}

void print_vect::print_itself(point pa, std::string str)
{
	for (int ii = 0; ii < pa.a.size(); ii++)
	{
		printf("a = %4d \t b = %4d \t x = %4.4f \t y = %4.4f \n", pa.a[ii], pa.b[ii], pa.x[ii], pa.y[ii]);
	}
	//printf("size: %d \t %s \n", pa.a.size(), str);
	std::cout << "size: " << pa.a.size() <<'\t' <<str << '\n';
}

void print_vect::print_itself(const std::vector<double> &x, const std::vector<double> &y, std::string str)
{
	for (int ii = 0; ii < x.size(); ii++)
	{
		printf("x1 = %4.4f \t y1 = %4.4f \n", x[ii], y[ii]);
	}
	std::cout << "size: " << x.size() << '\t' << str << '\n';
}