#ifndef getLogBins_h
#define getLogBins_h

#include <math.h>

void getLogBins(const float lower, const float higher, const int nBins, double bins[])
{
  float logBins[nBins+1];
  bins[0] = lower;
  bins[nBins] = higher;

  logBins[0] = std::log10(lower);
  logBins[nBins] = std::log10(higher);

  float interval = (logBins[nBins] - logBins[0])/nBins;

  for(int iter = 1; iter < nBins; iter++){
    logBins[iter] = logBins[0] + iter*interval;
    bins[iter] = std::pow(10., logBins[iter]);
  }

  return;
}

#endif
