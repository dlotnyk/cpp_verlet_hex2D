#include "coord.h"
#define PI 3.14159265359 // some useful constants
#define C120 cos(2*PI/3)
#define S120 sin(2*PI/3)

coord::coord()
{
	
	number_p = 1;
	std::vector<double> pp(number_p, 0);
	std::vector<double> ampl(number_p, 0);
	std::vector<double> orient(number_p, 0);
	beg.dt = 0;
	beg.eta = 0;
	beg.k = 0;
	beg.m = 0;
	beg.x0 = 0;
	beg.y0 = 0;
	//lattice_grid lat(number_p);+
	acceler.initialization(number_p);
	velo.initialization(number_p);
	std::cout << "default constructor1" << '\n';
}

coord::coord(const int &nn) : lattice_grid(nn)
{
	
	
	std::vector<double> pp(nn, 0);
	std::vector<double> ampl(nn, 0);
	std::vector<double> orient(nn, 0);
	number_p = nn;
	beg.dt = 0.1;
	beg.eta = 0;
	beg.k = 1;
	beg.m = 1;
	beg.x0 = 0.1;
	beg.y0 = 0.1;
	acceler.initialization(number_p);
	velo.initialization(number_p);
	//lattice_grid lat(number_p);
	//neighbor(1);
	change_lattice(beg.x0, beg.y0);
	std::cout << "default constructor2" << '\n';
}

coord::coord(const int &nn,const double &m, const double &k, const double &dt, const double &x0, const double &y0, const double &eta) : lattice_grid(nn)
{
	number_p = nn;
	std::vector<double> pp(number_p, 0);
	std::vector<double> ampl(number_p, 0);
	std::vector<double> orient(number_p, 0);
	beg.dt = dt;
	beg.eta = eta;
	beg.k = k;
	beg.m = m;
	beg.x0 = x0;
	beg.y0 = y0;
	acceler.initialization(number_p);
	velo.initialization(number_p);
	//lattice_grid lat(number_p);
	if (!beg.eta)
			change_lattice(); // only center vortex disported
		else
			change_lattice(beg.x0, beg.y0);
	
	std::cout << "default constructor3" << '\n';
	//lat.neighbor(0);
	
}

coord::coord(const int &nn, const double &m, const double &k, const double &dt, std::vector<double> &x0, std::vector<double> &y0, const double &eta) :lattice_grid(nn)
{
	number_p = nn;
	std::vector<double> pp(number_p, 0);
	std::vector<double> ampl(number_p, 0);
	std::vector<double> orient(number_p, 0);
	beg.dt = dt;
	beg.eta = eta;
	beg.k = k;
	beg.m = m;
	beg.x0 = 0;
	beg.y0 = 0;
	acceler.initialization(number_p);
	velo.initialization(number_p);
	//lattice_grid lat(number_p);
	change_lattice(x0, y0); // only center vortex disported
	std::cout << "default constructor4" << '\n';
}
coord::~coord()
{
	pp.clear();
	ampl.clear();
	orient.clear();
	acceler.clear();
	velo.clear();
	val.clear();
	//lattice.clear();
	//lat.~lattice_grid();
}

void coord::set_values(const std::vector<double> &xx, const std::vector<double> &yy) 
{
	pp.x = xx;
	pp.y = yy;
	print_itself(pp);
}

void coord::change_lattice(double x, double y)
{
	lattice.x[0] += x;
	lattice.y[0] += y;
	sample.x[0] = lattice.x[0];
	sample.y[0] = lattice.y[0];
}

void coord::change_lattice(double x, double y, int a)
{
	for (int ii = 0; ii < lattice.x.size(); ii++) {
		lattice.x[ii] += x;
		lattice.y[ii] += y;
		sample.x[ii] += x;
		sample.y[ii] += y;
	}
}

void coord::change_lattice()
{
	point dis;
	double randd;
	dis.initialization(lattice.a.size());
	for (int kk = 0; kk < dis.x.size(); kk++) {
		dis.x[kk] = rand() % 10;
		dis.y[kk] = rand() % 10;
		randd = sqrt(dis.x[kk] * dis.x[kk] + dis.y[kk] * dis.y[kk]);
		randd *= 10;
		dis.x[kk] /= randd;
		dis.y[kk] /= randd;
		lattice.x[kk] += dis.x[kk];
		lattice.y[kk] += dis.y[kk];
		sample.x[kk] += dis.x[kk];
		sample.y[kk] += dis.y[kk];
		//cout << "x= " << randn_x[kk] << '\t' << "y= " << randn_y[kk] << '\n';
	}

}

void coord::change_lattice(std::vector<double> &x, std::vector<double> &y)
{
	// x
	if (x.size() >= lattice.x.size())
	{
		for (int ii = 0; ii < lattice.x.size(); ii++)
		{
			lattice.x[ii] += x[ii];
			sample.x[ii] += x[ii];
		}
	}
	else
	{
		int s = 0, s1 = 0;
		while (s < lattice.x.size())
		{
			if (s1 < x.size())
			{
				lattice.x[s] += x[s1];
				sample.x[s] += x[s1];
				s1++;
			}
			if (s1 >= x.size())
				s1 = 0;
			s++;
		}
	}
	// y
	if (y.size() >= lattice.y.size())
	{
		for (int ii = 0; ii < lattice.y.size(); ii++)
		{
			lattice.y[ii] += y[ii];
			sample.y[ii] += y[ii];
		}
	}
	else
	{
		int s = 0, s1 = 0;
		while (s < lattice.y.size())
		{
			if (s1 < y.size())
			{
				lattice.y[s] += y[s1];
				sample.y[s] += y[s1];
				s1++;
			}
			if (s1 >= y.size())
				s1 = 0;
			s++;
		}
	}
}

void coord::change_velo()
{
	for (int ii = 0; ii < lattice.x.size(); ii++) {
		velo.x[ii] += 0.5*beg.dt*acceler.x[ii];
		velo.y[ii] += 0.5*beg.dt*acceler.y[ii];
	}
}

void coord::velo_begin(double vx, double vy)
{
	velo.x[0] += vx;
	velo.y[0] += vy;
}

void coord::velo_begin(double vx, double vy, int a)
{
	for (int ii = 0; ii < velo.x.size(); ii++) {
		velo.x[ii] += vx;
		velo.y[ii] += vy;
	}
}

void coord::change_coord() {
	double dt_s;
	dt_s = beg.dt*beg.dt;
	for (int ii = 0; ii < lattice.x.size(); ii++) {
		lattice.x[ii] = lattice.x[ii] + velo.x[ii] * beg.dt + 0.5*acceler.x[ii] * dt_s;
		lattice.y[ii] = lattice.y[ii] + velo.y[ii] * beg.dt + 0.5*acceler.y[ii] * dt_s;
		sample.x[ii] = sample.x[ii] + velo.x[ii] * beg.dt + 0.5*acceler.x[ii] * dt_s;
		sample.y[ii] = sample.y[ii] + velo.y[ii] * beg.dt + 0.5*acceler.y[ii] * dt_s;
		/*
		vx0[ii] = visco(vx0[ii]);
		vy0[ii] = visco(vy0[ii]);*/
	}
}

std::vector<double> coord::amplitude_point(const int &po) 
{
	printf("check size = %d \n", pp.x.size());
	if (pp.x.size() != pp.y.size() || po > pp.x.size())
	{
		printf("Size of x coordintaes and y should be equal \n");
		system("PAUSE");
		exit(EXIT_FAILURE);

	}
	ampl.clear();
	if (po != 0) {
		for (int n = 0; n < pp.x.size(); n++) {
			if ((po - 1) != n) {
				ampl.push_back(sqrt((pp.x[n] - pp.x[po - 1])*(pp.x[n] - pp.x[po - 1]) + (pp.y[n] - pp.y[po - 1])*(pp.y[n] - pp.y[po - 1])));
				printf("n= %d \n", n);
			}
		}
	}
	else {

		for (int n = 0; n < pp.x.size(); n++) {
			ampl.push_back(sqrt((pp.x[n])*(pp.x[n]) + (pp.y[n])*(pp.y[n])));
			printf("n1= %d \n", n);
		}
	}
	
	return ampl;
}

void coord::accelaration()
{
	//double xq[6], yq[6], xnq[6], ynq[6], xdis[6], ydis[6], xsum, ysum;
	point displace, unit_v, norm_disp;
	dot sum;
	displace.initialization(6);
	unit_v.initialization(6);
	norm_disp.initialization(6);
	
	//std::cout << "latt s = " << lattice.x.size() << '\n';

	for (int ff = 0; ff < lattice.x.size(); ff++) {
		neighbor(ff);
		// displacments
		for (int ii = 0; ii < 6; ii++) {
			displace.x[ii] = cent.x - neibor.x[ii];
			displace.y[ii] = cent.y - neibor.y[ii];

			//std::cout << "x[" << ii << "]=" << displace.x[ii] << "\t" << "y[" << ii << "]=" << displace.y[ii] << "\n";
		}
		// unit vectors
		for (int ii1 = 0; ii1 < 6; ii1++) {
			unit_v.x[ii1] = displace.x[ii1] / (sqrt(displace.x[ii1] * displace.x[ii1] + displace.y[ii1] * displace.y[ii1]));
			unit_v.y[ii1] = displace.y[ii1] / (sqrt(displace.x[ii1] * displace.x[ii1] + displace.y[ii1] * displace.y[ii1]));

			//std::cout << "xn[" << ii1 << "]=" << unit_v.x[ii1] << "\t" << "yn[" << ii1 << "]=" << unit_v.y[ii1] << "\n";
		}
		// norm disp
		for (int ii = 0; ii < 6; ii++) {
			norm_disp.x[ii] = displace.x[ii] - unit_v.x[ii];
			norm_disp.y[ii] = displace.y[ii] - unit_v.y[ii];

			//std::cout << "xdis[" << ii << "]=" << norm_disp.x[ii] << "\t" << "ydis[" << ii << "]=" << norm_disp.y[ii] << "\n";
		}
		// sum disp
		sum.x = 0;
		sum.y = 0;
		for (int ii = 0; ii < 6; ii++) {
			sum.x += norm_disp.x[ii];
			sum.y += norm_disp.y[ii];

		}
		//std::cout << "xsum= " << sum.x << "\t" << "ysum= " << sum.y << "\n";
		// acceleration
		acceler.x[ff] = -beg.k*sum.x / beg.m;
		acceler.y[ff] = -beg.k*sum.y / beg.m;

		//std::cout << "ax= " << acceler.x[ff] << "\t" << "ay= " << acceler.y[ff] << "\n";
	}
}

void coord::verlet()
{
	// velosity dt/2
	change_coord();
	change_velo();
	accelaration();
	change_velo();
	// velosity dt /2
}

void coord::get_values(int num)
{
	val.initialization(num);
	for (int ii = 0; ii < val.x.size(); ii++) {
		val.x[ii] = sample.x[ii];
		val.y[ii] = sample.y[ii];
		val.vx[ii] = velo.x[ii];
		val.vy[ii] = velo.y[ii];
		val.ax[ii] = acceler.x[ii];
		val.ay[ii] = acceler.y[ii];
	}
}
#undef S120
#undef C120
#undef PI