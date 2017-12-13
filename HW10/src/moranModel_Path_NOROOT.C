#include <iostream>
#include <vector>
#include <fstream>
#include <ctime>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real_distribution.hpp>

#include "include/doGlobalDebug.h"
#include "include/getLinBins.h"
#include "include/getLogBins.h"
#include "include/cppWatch.h"

std::string prettyString(const double inVal, const int prec, const bool doDot)
{
  std::string retStr = std::to_string(inVal);
  while(retStr.find(".") < retStr.size()-1-prec){retStr.replace(retStr.size()-1, 1,"");}
  if(doDot) retStr.replace(retStr.find("."), 1, "p");
  return retStr;
}

int moranModel_Path_NOROOT(const double s0A, const double s0B, const double s0AB, const unsigned int nSim = 500, const int nPopSize=10, const double overrideMuteRate = -1)
{
  boost::random::mt19937 gen(std::time(0));
  boost::random::uniform_real_distribution<> dist(0, 1);
  time_t now = time(0);
  tm *ltm = localtime(&now);

  std::string fullDateStr = std::to_string(ltm->tm_year + 1900);
  if(ltm->tm_mon + 1 < 10) fullDateStr = fullDateStr + "0";
  fullDateStr = fullDateStr + std::to_string(ltm->tm_mon + 1);
  if(ltm->tm_mday < 10) fullDateStr = fullDateStr + "0";
  fullDateStr = fullDateStr + std::to_string(ltm->tm_mday);

  //decalare mutation rates to scan, to reduce/increase number of points, tweak 'nMuteRate'
  const double mutRateLow = 0.0001;
  const double mutRateHi = 1.;
  
  int tempNMuteRate = 20;
  if(overrideMuteRate > 0) tempNMuteRate = 1;
  const int nMuteRate = tempNMuteRate;
  double mutRate[nMuteRate+1];
  mutRate[0] = overrideMuteRate;
  if(overrideMuteRate < 0) getLogBins(mutRateLow, mutRateHi, nMuteRate, mutRate);

  std::string csvFileName = "output/histFile_ProbAandB_sA" + prettyString(s0A, 3, true) + "_sB" + prettyString(s0B, 3, true) + "_sAB" + prettyString(s0AB, 3, true) + "_nSim" + std::to_string(nSim) + "_nPopSize" + std::to_string(nPopSize) + "_OverrideMuteRate" + prettyString(overrideMuteRate, 6, true) + "_" + fullDateStr  + ".csv";

  std::ofstream fileA(csvFileName.c_str());
  fileA << "muteRate(sAB=" << s0AB << "),Prob. Path A(sA=" << s0A << "),Error on Prob. Path A,Prob. Path B(sB=" << s0B << "),Error on Prob. Path B,";
  fileA << "StartA/FixA,Error StartA/FixA,StartA/FixAB,Error StartA/FixAB,StartA/FixB,Error StartA/FixB,StartA/FixBA,Error StartA/FixBA,";
  fileA << "StartB/FixA,Error StartB/FixA,StartB/FixAB,Error StartB/FixAB,StartB/FixB,Error StartB/FixB,StartB/FixBA,Error StartB/FixBA,";
  fileA << std::endl;
  fileA.close();


  for(int muI = 0; muI < nMuteRate+1; ++muI){
    if(overrideMuteRate > 0 && muI != 0) break;

    std::cout << "Running mutation " << muI << "/" << nMuteRate << ": " << mutRate[muI] << std::endl;

    std::vector<std::string> sGen_Win_;
    std::vector<std::string> sGen_WinMethod_;
    sGen_Win_.reserve(nSim);
    sGen_WinMethod_.reserve(nSim);
    int counter = 0;

    cppWatch extrap;
    extrap.start();
    bool timerIsDone = false;
    const unsigned int winnersToEstimate = 20;
    unsigned int prevWinners = 0;
   
    while(nSim > sGen_Win_.size()){
      if(sGen_Win_.size() >= winnersToEstimate && !timerIsDone){
	extrap.stop();
	std::cout << " Time for first " << winnersToEstimate << "  winners: " << 1000.0 * extrap.total() / CLOCKS_PER_SEC << " ms" << std::endl;
	std::cout << " First " << winnersToEstimate << " winners: ";
	for(unsigned int i = 0; i < sGen_Win_.size(); ++i){
	  std::cout << sGen_Win_.at(i) << ",";
	}
	std::cout << std::endl;
	std::cout << " Extrap. time left: " << (1000.0 * extrap.total() / CLOCKS_PER_SEC)*(double(nSim - winnersToEstimate))/((double)winnersToEstimate) << " ms" << std::endl;;
	timerIsDone = true;
	extrap.start();
      }

      if(timerIsDone && sGen_Win_.size()%winnersToEstimate == 0 && sGen_Win_.size() != prevWinners){
	std::cout << " Now at " << sGen_Win_.size() << "/" << nSim << " winners" << std::endl;
	prevWinners = sGen_Win_.size();
      }

      std::vector<int> pop;
      std::vector<double> relFit;     

      pop.push_back(nPopSize-1);
      pop.push_back(0);
      pop.push_back(0);
      pop.push_back(0);
      pop.push_back(0);

      if(counter%2 == 0) pop.at(1) = 1;
      else pop.at(2) = 1;
      
      relFit.push_back(1.);
      relFit.push_back(1. + s0A);
      relFit.push_back(1. + s0B);
      relFit.push_back(1. + s0AB);
      relFit.push_back(1. + s0AB);

      int pops = 1;

      while(pops != nPopSize){
	double denom = 0.;
	for(unsigned int jI = 0; jI < pop.size(); ++jI){
	  denom += pop.at(jI)*relFit.at(jI);
	}
	
	std::vector<double> probRep;
	std::vector<double> probDie;

	for(unsigned int jI = 0; jI < pop.size(); ++jI){
	  probRep.push_back(pop.at(jI)*relFit.at(jI)/denom);
	  probDie.push_back(double(pop.at(jI))/(double)nPopSize);
	}
	
	const int nMatrix = probRep.size();
	double matrix[nMatrix][nMatrix];
	double newDenom = 0.;
	
	for(int mI = 0; mI < nMatrix; ++mI){
	  for(int nI  = 0; nI < nMatrix; ++nI){
	    matrix[mI][nI] = probRep.at(mI)*probDie.at(nI);
	    if(mI != nI) newDenom += matrix[mI][nI];
	  }
	}

	for(int mI = 0; mI < nMatrix; ++mI){
	  for(int nI  = 0; nI < nMatrix; ++nI){
	    matrix[mI][nI] /= newDenom;
	  }
	}

	double draw1 = dist(gen);
	
	double fullProb = 0;
	int repPos = -1;
	int deathPos = -1;
	bool isFound = false;
	
	for(int mI = 0; mI < nMatrix; ++mI){
	  for(int nI  = 0; nI < nMatrix; ++nI){
	    if(mI == nI) continue;
	    
	    fullProb += matrix[mI][nI];
	    if(draw1 < fullProb){
	      repPos = mI;
	      deathPos = nI;
	      isFound = true;
	      break;
	    }
	  }
	  if(isFound) break;
	}

	pop.at(deathPos) -= 1;
	double draw2 = dist(gen);
      
	if(draw2 < mutRate[muI] && repPos != 3 && repPos != 4){
	  if(repPos == 1) pop.at(3) += 1;
	  else if(repPos == 2) pop.at(4) += 1;
	  else{
	    double draw3 = dist(gen);
	    if(draw3 < .5) pop.at(2) += 1;
	    else pop.at(1) += 1;
	  }
	}
	else pop.at(repPos) += 1;

	bool doDraw = pop.at(0) == 0 && pop.at(1) == 0 && pop.at(2) == 0 && pop.at(3) != 0 && pop.at(4) != 0;	
	if(doDraw){
	  double draw4 = dist(gen);

	  if(draw4 < ((double)pop.at(3))/((double)nPopSize)){
	    pop.at(3) = nPopSize;
	    pop.at(4) = 0;
	  }
	  else{
	    pop.at(3) = 0;
	    pop.at(4) = nPopSize;
	  }
	}

	pops = pop.at(0);
	if(pops < pop.at(1) + pop.at(3)) pops = pop.at(1) + pop.at(3);
	else if(pops < pop.at(2) + pop.at(4)) pops = pop.at(2) + pop.at(4);
      }
    
      std::string winPathStr = "b";
      if(counter%2 == 0) winPathStr = "a";
      ++counter;
    
      if(pop.at(0) == nPopSize) continue;

      if(pop.at(1) + pop.at(3) == nPopSize){
	sGen_Win_.push_back("a");

	if(pop.at(1) == nPopSize) winPathStr = winPathStr + "_a";
	else winPathStr = winPathStr + "_ab";

	sGen_WinMethod_.push_back(winPathStr);
      }
      else if(pop.at(2) + pop.at(4) == nPopSize){
	sGen_Win_.push_back("b");

	if(pop.at(2) == nPopSize) winPathStr = winPathStr + "_b";
	else winPathStr = winPathStr + "_ba";

	sGen_WinMethod_.push_back(winPathStr);
      }
    }

    extrap.stop();

    std::cout << " Time for total " << nSim << "  winners: " << 1000.0 * extrap.total() / CLOCKS_PER_SEC << " ms" << std::endl;

    
    double aCounts = 0;
    double bCounts = 0;

    double aaCounts = 0;
    double aabCounts = 0;
    double abCounts = 0;
    double abaCounts = 0;

    double baCounts = 0;
    double babCounts = 0;
    double bbCounts = 0;
    double bbaCounts = 0;

    for(unsigned int i = 0; i < sGen_Win_.size(); ++i){
      if(sGen_Win_.at(i).find("a") != std::string::npos && sGen_Win_.at(i).size() == 1) ++aCounts;
      else if(sGen_Win_.at(i).find("b") != std::string::npos && sGen_Win_.at(i).size() == 1) ++bCounts;

      if(sGen_WinMethod_.at(i).find("a_a") != std::string::npos && sGen_WinMethod_.at(i).size() == 3) ++aaCounts;
      else if(sGen_WinMethod_.at(i).find("a_ab") != std::string::npos && sGen_WinMethod_.at(i).size() == 4) ++aabCounts;
      else if(sGen_WinMethod_.at(i).find("a_b") != std::string::npos && sGen_WinMethod_.at(i).size() == 3) ++abCounts;
      else if(sGen_WinMethod_.at(i).find("a_ba") != std::string::npos && sGen_WinMethod_.at(i).size() == 4) ++abaCounts;
      else if(sGen_WinMethod_.at(i).find("b_a") != std::string::npos && sGen_WinMethod_.at(i).size() == 3) ++baCounts;
      else if(sGen_WinMethod_.at(i).find("b_ab") != std::string::npos && sGen_WinMethod_.at(i).size() == 4) ++babCounts;
      else if(sGen_WinMethod_.at(i).find("b_b") != std::string::npos && sGen_WinMethod_.at(i).size() == 3) ++bbCounts;
      else if(sGen_WinMethod_.at(i).find("b_ba") != std::string::npos && sGen_WinMethod_.at(i).size() == 4) ++bbaCounts;
    }
    
    fileA.open(csvFileName.c_str(), std::ios_base::app);
    fileA << mutRate[muI] << "," << aCounts/(double)nSim << "," << std::sqrt(aCounts)/(double)nSim << "," << bCounts/(double)nSim << "," << std::sqrt(bCounts)/(double)nSim << ",";
    fileA << aaCounts/(double)nSim << "," << std::sqrt(aaCounts)/(double)nSim << ",";
    fileA << aabCounts/(double)nSim << "," << std::sqrt(aabCounts)/(double)nSim << ",";
    fileA << abCounts/(double)nSim << "," << std::sqrt(abCounts)/(double)nSim << ",";
    fileA << abaCounts/(double)nSim << "," << std::sqrt(abaCounts)/(double)nSim << ",";
    fileA << baCounts/(double)nSim << "," << std::sqrt(baCounts)/(double)nSim << ",";
    fileA << babCounts/(double)nSim << "," << std::sqrt(babCounts)/(double)nSim << ",";
    fileA << bbCounts/(double)nSim << "," << std::sqrt(bbCounts)/(double)nSim << ",";
    fileA << bbaCounts/(double)nSim << "," << std::sqrt(bbaCounts)/(double)nSim << ",";
    fileA << std::endl;

    fileA.close();
  }

  return 0;
}


int main(int argc, char* argv[])
{
  if(argc != 4 && argc != 5 && argc != 6 && argc != 7){
    std::cout << "Usage: ./moranModel_Path_NOROOT.exe <s0A> <s0B> <s0AB> <nSim-optional> <nPopSize-optional>" << std::endl;
  }

  int retVal = 0;
  if(argc == 4) retVal += moranModel_Path_NOROOT(std::stof(argv[1]), std::stof(argv[2]), std::stof(argv[3]));
  else if(argc == 5) retVal += moranModel_Path_NOROOT(std::stof(argv[1]), std::stof(argv[2]), std::stof(argv[3]), std::stoi(argv[4]));
  else if(argc == 6) retVal += moranModel_Path_NOROOT(std::stof(argv[1]), std::stof(argv[2]), std::stof(argv[3]), std::stoi(argv[4]), std::stoi(argv[5]));
  else if(argc == 7) retVal += moranModel_Path_NOROOT(std::stof(argv[1]), std::stof(argv[2]), std::stof(argv[3]), std::stoi(argv[4]), std::stoi(argv[5]), std::stof(argv[6]));
  return retVal;
}
