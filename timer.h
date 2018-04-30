#ifndef TIMER_H_
#define TIMER_H_
#include <chrono>
class timer
{
public:
	double time;
	timer();
	~timer();
	double getClock();
};

#endif TIMER_H_