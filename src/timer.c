#include <time.h>

#include "timer.h"

double elapsed_time(clock_t s_time, clock_t e_time)
{
	return ((double) (e_time - s_time))/((double) CLOCKS_PER_SEC);
}
