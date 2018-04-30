#include "coord.h"
#define NUMBER_OF_STEPS 131072
//#include "lattice_grid.h"
int main() {

	//initialization
	timer t1; // starting timer
	int number = 19; // number of oscillationg vorticies
	int save_num = 7; // number of them to save
	double m = 1;
	double k = 10;
	double dt = 0.01;
	double eta=0;
	double curr_time = 0;
	std::vector<double> x = { 0.1, 0.15, 0 }; // displacements
	std::vector<double> y = { 0.1, 0.3, 0 };
	print_vect CX;

	coord BX(number, m, k, dt, 0.1, 0.1, eta); // first timestep
	BX.accelaration();

	FILE *file; // file
	int err;
	//filename = string("disp_s%u_num%u.dat",s,vva);		
	err = fopen_s(&file, "disp.dat", "w");
	if (err == 0)
	{
		printf("The file 'disp.dat' was opened\n");
		for (int jj = 0; jj < NUMBER_OF_STEPS+1; jj++) {
			curr_time = jj*dt;
			BX.get_values(save_num);
			for (int ii = 0; ii < save_num; ii++) {
				fprintf(file, "%.6f \t %.6f \t %.6f \t %.6f \t %.6f \t %.6f \t %.6f", curr_time, BX.val.x[ii], BX.val.y[ii], BX.val.vx[ii], BX.val.vy[ii], BX.val.ax[ii], BX.val.ay[ii]);
				fprintf(file, "\n");
			}
			BX.verlet();
		}
	}
	else
	{
		printf("The file 'disp.dat' was not opened\n");
	}
	fclose(file);
	

	timer t2; // stop the timer
	std::cout << "operation time = " << t2.time - t1.time << " sec \n";
	system("PAUSE");
	return 0;
};