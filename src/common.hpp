#define INVLOG2 1.44269504089
#include<boost/filesystem.hpp>
#include <my_sort.hpp>
#include <time.h>

using namespace std;
#define MAXCLUSTER 10000
#define NCOLORS 500
#define NCOLORS1OFF 501
#define MINBEADS 20
#define FLOATNULL std::numeric_limits<float>::quiet_NaN()
#define DOUBLENULL std::numeric_limits<double>::quiet_NaN()
#define MIXTURE_RATIO .667
#define MAX_MIXTURE_RATIO .700
#define MIXTURE_STD .015
//timers

struct timespec start_time;
struct timespec end_time;

double get_elapsed_time(const struct timespec *start_time, const struct timespec *end_time)
{
    int64_t sec = end_time->tv_sec - start_time->tv_sec;
    int64_t nsec;
    if (end_time->tv_nsec >= start_time->tv_nsec) {
        nsec = end_time->tv_nsec - start_time->tv_nsec;
    } else {
        nsec = 1000000000 - (start_time->tv_nsec - end_time->tv_nsec);
        sec -= 1;
    }
    return ((double) sec + (double) nsec *1e-9);  
}
