/*
  Developed by Sandeep Sharma 
  Copyright (c) 2017, Sandeep Sharma
  
  This file is part of DICE.
  
  This program is free software: you can redistribute it and/or modify it under the terms
  of the GNU General Public License as published by the Free Software Foundation, 
  either version 3 of the License, or (at your option) any later version.
  
  This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
  without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  
  See the GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License along with this program. 
  If not, see <http://www.gnu.org/licenses/>.
*/
#include <algorithm>
#include <random>
#include <chrono>
#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <boost/format.hpp>
#include <boost/algorithm/string.hpp>
#ifndef SERIAL
//#include "mpi.h"
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi.hpp>
#endif
#include "evaluatePT.h"
#include "Determinants.h"
#include "CPSSlater.h"
#include "HFWalker.h"
#include "AGP.h"
//#include "AGPWalker.h"
#include "input.h"
#include "integral.h"
#include "SHCIshm.h"
#include "math.h"
#include "Profile.h"
#include "CIWavefunction.h"
#include "runVMC.h"
#include "Lanczos.h"

using namespace Eigen;
using namespace boost;
using namespace std;


int main(int argc, char* argv[]) {

#ifndef SERIAL
  boost::mpi::environment env(argc, argv);
  boost::mpi::communicator world;
#endif
  startofCalc = getTime();

  initSHM();
  //license();

  string inputFile = "input.dat";
  if (argc > 1)
    inputFile = string(argv[1]);
  readInput(inputFile, schd);

  generator = std::mt19937(schd.seed+commrank);

  readIntegralsAndInitializeDeterminantStaticVariables("FCIDUMP");

  double ham, Ept2, stddev, stddev2, rk;
  if (schd.wavefunctionType == "CPSSlater") {
    CPSSlater<CPS, Slater> wave; HFWalker<CPS, Slater> walk;
    wave.readWave(); wave.initWalker(walk);
    getStochasticEnergyContinuousTime(wave, walk, ham, stddev, rk, schd.stochasticIter, 1.e-5);
    Ept2 = evaluatePTStochastic(wave, walk, schd.PTlambda, ham,
                                     stddev2, rk,
                                     schd.stochasticIter, 0.5e-3);
  }
  else if (schd.wavefunctionType == "JastrowSlater") {
    CPSSlater<Jastrow, Slater> wave; HFWalker<Jastrow, Slater> walk;
    wave.readWave();
    wave.initWalker(walk);
    getStochasticEnergyContinuousTime(wave, walk, ham, stddev, rk, schd.stochasticIter, 1.e-5);
    if (commrank == 0) cout << ham <<endl;
    Ept2 = evaluatePTDeterministic(wave, walk, ham,
                                   stddev2, rk,
                                   schd.stochasticIter, 0.5e-3);
    cout << Ept2<<endl;
      Ept2 = evaluatePTStochastic(wave, walk, schd.PTlambda, ham,
                                  stddev2, rk,
                                     schd.stochasticIter, 0.5e-3);
    if (commrank == 0) cout << Ept2 <<endl;
  }
  else if (schd.wavefunctionType == "LanczosJastrowSlater") {
    Lanczos<CPSSlater<Jastrow, Slater>> wave; HFWalker<Jastrow, Slater> walk;
    wave.readWave();
    wave.initWalker(walk);
    getStochasticEnergyContinuousTime(wave, walk, ham, stddev, rk, schd.stochasticIter, 1.e-5);
    if (commrank == 0) cout << ham <<endl;
    if (schd.deterministic)
      Ept2 = evaluatePTDeterministic(wave, walk, ham,
                                       stddev2, rk,
                                       schd.stochasticIter, 0.5e-3);
    else
      Ept2 = evaluatePTStochastic(wave, walk, schd.PTlambda, ham,
                                       stddev2, rk,
                                       schd.stochasticIter, 0.5e-3);

    if (commrank == 0) cout << Ept2 <<endl;
  }

  if (commrank == 0) {
    std::cout << format("%14.8f (%8.2e) %14.8f (%8.2e) %10.2f %10i %8.2f\n") 
      %ham % stddev %(ham+Ept2) %stddev2 %rk %(schd.stochasticIter) %( (getTime()-startofCalc));

  }


  boost::interprocess::shared_memory_object::remove(shciint2.c_str());
  boost::interprocess::shared_memory_object::remove(shciint2shm.c_str());
  return 0;
}
