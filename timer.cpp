#include "timer.h"



timer::timer()
{
	time = getClock();
}

timer::~timer()
{
}
double timer::getClock()
{
		auto now = std::chrono::high_resolution_clock::now();
		return std::chrono::duration_cast<std::chrono::microseconds>(now.time_since_epoch()).count() * 1e-6;
	
}