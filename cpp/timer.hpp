/*
 * Thread safe, high resolution monotonic timing routines.
 *
 * You to link with '-lrt'.
 */


#ifndef _PALEON_TIMER_HPP_
#define _PALEON_TIMER_HPP_

#include <cfloat>
#include <iostream>
#include <vector>

#include <time.h>

using namespace std;

namespace paleon
{
  class Timer
  {
    struct timer {
      struct timespec start, end;
    };

    vector<struct timer> timers;
    vector<double> seconds;

  public:
    Timer(int nthreads) : timers(nthreads), seconds(nthreads) { }

    void tic(int thread)
    {
      clock_gettime(CLOCK_MONOTONIC, &timers.at(thread).start);
    }

    void toc(int thread)
    {
      clock_gettime(CLOCK_MONOTONIC, &timers.at(thread).end);

      struct timespec _tic = timers.at(thread).start;
      struct timespec _toc = timers.at(thread).end;
      seconds.at(thread) = (_toc.tv_sec - _tic.tv_sec) + (_toc.tv_nsec - _tic.tv_nsec)*1e-9;
    }

    void echo(string msg)
    {
      double min = DBL_MAX;
      double max = DBL_MIN;
      double tot = 0.0;

      for (int i=0; i<seconds.size(); i++) {
	double x = seconds[i];
	if (x < min) min = x;
	if (x > max) max = x;
	tot += x;
      }
      double avg = tot / seconds.size();

      cout << msg << " min/max/avg/tot: " << min << " " << max << " " << avg << " " << tot << endl;
      // XXX: span
    }

  };
}


#endif
