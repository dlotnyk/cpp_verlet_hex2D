#include "lattice_grid.h"
#define PI 3.14159265359 // some useful constants
#define C120 cos(2*PI/3)
#define S120 sin(2*PI/3)

lattice_grid::lattice_grid()
{
	number_p = 1;
	base_hex_lattice.initialization(6);
	
	base_hex_lattice.a = { 0, 1, 1, 0, -1, -1 };
	base_hex_lattice.b = { 1, 1, 0, -1, -1, 0 };

	lattice.initialization(6);
	neibor.initialization(6);
	

	hex_transform(base_hex_lattice);

	gen_lattice();
	gen_sample();

	std::cout << "default constructor0" << '\n';
}

lattice_grid::lattice_grid(const int &nn)
{
	number_p = nn;
	base_hex_lattice.initialization(6);

	base_hex_lattice.a = { 0, 1, 1, 0, -1, -1 };
	base_hex_lattice.b = { 1, 1, 0, -1, -1, 0 };

	lattice.initialization(number_p);
	neibor.initialization(6);

	hex_transform(base_hex_lattice);
	//print_itself(base_hex_lattice, "base vortecies");

	gen_lattice();
	gen_sample();
	//neighbor(1);
}

lattice_grid::~lattice_grid()
{
	base_hex_lattice.clear();
	lattice.clear();
	neibor.clear();
	sample.clear();
}

// hexahonal to cartesian transfordm
void lattice_grid::hex_transform(point &hextocart)
{
	if (hextocart.a.size() != hextocart.x.size())
	{
		hextocart.x = std::vector<double>(hextocart.a.size(), 0);
		hextocart.y = std::vector<double>(hextocart.a.size(), 0);
	}
	for (int ii = 0; ii < hextocart.a.size(); ii++)
	{
		hextocart.x[ii] = hextocart.a[ii] + hextocart.b[ii] * C120;
		hextocart.y[ii] = hextocart.b[ii] * S120;
	}
	//print_itself(hextocart);
}

// generate vibrating points
void lattice_grid::gen_lattice() // only vibrating
{
	// generate a points in hex lattice base_hex_lattice
	point trash;
	int num_count = 1;
	int hex_count = 1;
	//int trash_count = 0;
	lattice.a[0] = 0;
	lattice.b[0] = 0;
	for (int kk = 1; kk < number_p; kk++)
	{

		while (!trash.a.empty())
		{
			//printf(" k = %d trash = %d \n", kk, trash.a.size());
			if (num_count == number_p) break;
			lattice.a[num_count] = trash.a[0];
			lattice.b[num_count] = trash.b[0];
			num_count++;
			trash.a.erase(trash.a.begin());
			trash.b.erase(trash.b.begin());
		}
		if (num_count == number_p) break;

		if (num_count == number_p) break;

		for (int ii = -kk; ii < kk + 1; ii++)
		{
			//printf("i = %d \t k = %d \n", ii, kk);
			if (num_count == number_p) break;
			for (int jj = -kk; jj < kk + 1; jj++)
			{
				//printf("j = %d \n", jj);
				if (num_count == number_p) break;
				if (abs(ii) == kk || abs(jj) == kk)
				{
					if (ii*jj >= 0)
					{
						lattice.a[num_count] = ii;
						lattice.b[num_count] = jj;
						num_count++;
					}
					else
					{
						trash.a.push_back(ii);
						trash.b.push_back(jj);

					}
				}
			}

		}
	}
	std::string str1("lattice");
	hex_transform(lattice);
	//print_itself(lattice, str1);
	std::string str2("trash");
	hex_transform(trash);
	//print_itself(trash, str2);
	trash.clear();
}

// generate all points
void lattice_grid::gen_sample() {

	point an;

	for (int ii = 0; ii < number_p; ii++) {
		for (int jj = 0; jj < 6; jj++) {
			an.a.push_back(lattice.a[ii] + base_hex_lattice.a[jj]);
			an.b.push_back(lattice.b[ii] + base_hex_lattice.b[jj]);
		}
	}
	hex_transform(an);
	sample.a = lattice.a;
	sample.b = lattice.b;
	for (int ii = 0; ii < an.a.size(); ii++) {
		bool coinc(0);
		for (int jj = 0; jj < sample.a.size(); jj++) {
			if (an.a[ii] == sample.a[jj] && an.b[ii] == sample.b[jj])
				coinc = 1;
		}
		if (!coinc) {
			sample.a.push_back(an.a[ii]);
			sample.b.push_back(an.b[ii]);
		}
	}
	sample.x = std::vector<double>(sample.a.size(), 0);
	sample.y = std::vector<double>(sample.a.size(), 0);

	hex_transform(sample);
	//print_itself(sample, "all points");
	an.clear();
}

// neighbours to numbered point
void lattice_grid::neighbor(int ini_point)
{
	cent.x = lattice.x[ini_point];
	cent.y = lattice.y[ini_point];
	for (int jj = 0; jj < 6; jj++) {
		for (int ii = 0; ii < sample.a.size(); ii++)
		{
			if (sample.a[ii] == (base_hex_lattice.a[jj] + lattice.a[ini_point]) && sample.b[ii] == (base_hex_lattice.b[jj] + lattice.b[ini_point]))
			{
				neibor.a[jj] = sample.a[ii];
				neibor.b[jj] = sample.b[ii];
				neibor.x[jj] = sample.x[ii];
				neibor.y[jj] = sample.y[ii];
				//cout << "x0 = " << xq0[jj] <<'\t'<< "y0= "<< yq0[jj]<<'\n';
			}
		}
	}
	//hex_transform(neibor);
	//print_itself(neibor, "neighbours");
}


#undef S120
#undef C120
#undef PI