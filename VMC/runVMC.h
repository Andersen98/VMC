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

#pragma once

#include <sys/stat.h>
#include "input.h"
#include "evaluateE.h"
#include "amsgrad.h"
#include "sgd.h"
#include "sr.h"
#include "linearMethod.h"

using functor1 = boost::function<void (VectorXd&, VectorXd&, double&, double&, double&)>;
using functor2 = boost::function<void (VectorXd&, VectorXd&, VectorXd&, DirectMetric&, double&, double&, double&)>;
using functor3 = boost::function<void (VectorXd&, VectorXd&, double&, double&, double&)>;
using functor4 = boost::function<void (VectorXd&, VectorXd&, VectorXd&, DirectMetric&, double&, double&, double&)>;
using functor5 = boost::function<void (VectorXd&, VectorXd&, MatrixXd&, MatrixXd&, double&, double&, double&)>;


template<typename Wave, typename Walker>
void runVMC(Wave& wave, Walker& walk) {

  if (schd.restart) wave.readWave();
  VectorXd vars; wave.getVariables(vars);
  getGradientWrapper<Wave, Walker> wrapper(wave, walk, schd.stochasticIter, schd.ctmc);
  functor1 getStochasticGradient = boost::bind(&getGradientWrapper<Wave, Walker>::getGradient, &wrapper, _1, _2, _3, _4, _5, schd.deterministic);
  functor2 getStochasticGradientMetric = boost::bind(&getGradientWrapper<Wave, Walker>::getMetric, &wrapper, _1, _2, _3, _4, _5, _6, _7, schd.deterministic);

  if (schd.method == amsgrad || schd.method == amsgrad_sgd) {
    AMSGrad optimizer(schd.stepsize, schd.decay1, schd.decay2, schd.maxIter, schd.avgIter);
    optimizer.optimize(vars, getStochasticGradient, schd.restart);
  }
  else if (schd.method == sgd) {
    SGD optimizer(schd.stepsize, schd.maxIter);
    optimizer.optimize(vars, getStochasticGradient, schd.restart);
  }
  else if (schd.method == sr) {
    SR optimizer(schd.stepsize, schd.maxIter);
    optimizer.optimize(vars, getStochasticGradientMetric, schd.restart);
  }
  else if (schd.method == linearmethod) {
    
  }
  
}


template<typename Wave, typename Walker>
void runVMCRealSpace(Wave& wave, Walker& walk) {

  if (schd.restart) wave.readWave();
  VectorXd vars; wave.getVariables(vars);
  getGradientWrapper<Wave, Walker> wrapper(wave, walk, schd.stochasticIter, schd.ctmc);

  functor3 getStochasticGradientRealSpace = boost::bind(&getGradientWrapper<Wave, Walker>::getGradientRealSpace, &wrapper, _1, _2, _3, _4, _5, schd.deterministic);
  functor4 getStochasticGradientMetricRealSpace = boost::bind(&getGradientWrapper<Wave, Walker>::getMetricRealSpace, &wrapper, _1, _2, _3, _4, _5, _6, _7, schd.deterministic);
  functor5 getStochasticGradientHessianRealSpace = boost::bind(&getGradientWrapper<Wave, Walker>::getHessianRealSpace, &wrapper, _1, _2, _3, _4, _5, _6, _7, schd.deterministic);

  if (schd.walkerBasis == REALSPACESTO || schd.walkerBasis == REALSPACEGTO) {
    if (schd.method == amsgrad || schd.method == amsgrad_sgd) {
      AMSGrad optimizer(schd.stepsize, schd.decay1, schd.decay2, schd.maxIter, schd.avgIter);
      optimizer.optimize(vars, getStochasticGradientRealSpace, schd.restart);
    }
    else if (schd.method == sgd) {
      SGD optimizer(schd.stepsize, schd.maxIter);
      optimizer.optimize(vars, getStochasticGradientRealSpace, schd.restart);
    }
    else if (schd.method == sr) {
      SR optimizer(schd.stepsize, schd.maxIter);
      optimizer.optimize(vars, getStochasticGradientMetricRealSpace, schd.restart);
    }
    else if (schd.method == linearmethod) {
      LM optimizer(schd.stepsize, schd.maxIter);
      optimizer.optimize(vars, getStochasticGradientHessianRealSpace, schd.restart);
    }
  }
}

