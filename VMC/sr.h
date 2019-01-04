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
#ifndef OPTIMIZERSR_HEADER_H
#define OPTIMIZERSR_HEADER_H
#include <Eigen/Dense>
#include <boost/serialization/serialization.hpp>
#include "iowrapper.h"
#include "global.h"
#include <boost/format.hpp>
#include <boost/algorithm/string.hpp>
#include <fstream>

#ifndef SERIAL
#include "mpi.h"
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi.hpp>
#endif

using namespace Eigen;
using namespace boost;
using namespace std;

class DirectMetric
{
    public:
    VectorXd T;
    MatrixXd Vectors;
    //std::vector<double> T;
    //std::vector<VectorXd> Vectors;
    //VectorXd H;
    double diagshift;
    MatrixXd S;
    DirectMetric(double _diagshift) : diagshift(_diagshift) {}

    void BuildMetric()
    {
      double Tau = 0.0;
      int dim = Vectors.col(0).rows();
      S = MatrixXd::Zero(dim, dim);
      for (int i = 0; i < Vectors.cols(); i++)
      {
        S += T(i) * Vectors.col(i) * Vectors.col(i).adjoint();
        Tau += T(i);
      }

#ifndef SERIAL
      MPI_Allreduce(MPI_IN_PLACE, &(Tau), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, &(S(0,0)), S.rows() * S.cols(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif

      S /= Tau; 
/*
      if (commrank == 0) std::cout << S << std::endl << std::endl;
      if (commrank == 1) std::cout << S << std::endl << std::endl;
*/
    }
      
    
  void multiply(VectorXd &x, VectorXd& Ax)
  {
    double Tau = 0.0;
    int dim = x.rows();
    if (Ax.rows() != x.rows())
      Ax = VectorXd::Zero(dim);
    else
      Ax.setZero();
    for (int i = 0; i < Vectors.cols(); i++)
    {
      //double factor = x.dot(Vectors[i]) * T[i];
      double factor = Vectors.col(i).adjoint() * x;
      Ax += T(i) * Vectors.col(i) * factor;
      Tau += T[i];
    }

#ifndef SERIAL
      MPI_Allreduce(MPI_IN_PLACE, &(Ax(0)), Ax.rows(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, &(Tau), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
      Ax /= Tau;

      for (int i = 1; i < Ax.rows(); i++)
        Ax(i) += diagshift * x(i);
  }
};

void ConjGrad(DirectMetric &A, VectorXd &b, int n, VectorXd &x);

void PInv(MatrixXd &A, MatrixXd &Ainv);

class SR
{
  private:
    friend class boost::serialization::access;
    template <class Archive>
    void serialize(Archive &ar, const unsigned int version)
    {
        ar & stepsize
            & iter;
    }

  public:
    double stepsize;
    int maxIter;
    int iter;

    SR(double pstepsize=0.001, int pmaxIter=1000) : stepsize(pstepsize), maxIter(pmaxIter)
    {
        iter = 0;
    }

    void write(VectorXd& vars)
    {
        if (commrank == 0)
        {
	    char file[5000];
            sprintf(file, "sr.bkp");
            std::ofstream ofs(file, std::ios::binary);
            boost::archive::binary_oarchive save(ofs);
            save << *this;
            save << vars;
            ofs.close();
        }
    }

    void read(VectorXd& vars)
    {
        if (commrank == 0)
        {
	    char file[5000];
            sprintf(file, "sr.bkp");
            std::ifstream ifs(file, std::ios::binary);
	    boost::archive::binary_iarchive load(ifs);
            load >> *this;
            load >> vars;
            ifs.close();
        }
    }

   template<typename Function>
   void optimize(VectorXd &vars, Function& getMetric, bool restart)
   {
     if (restart)
     {
       if (commrank == 0) 
         read(vars);
#ifndef SERIAL
	    boost::mpi::communicator world;
	    boost::mpi::broadcast(world, *this, 0);
	    boost::mpi::broadcast(world, vars, 0);
#endif
     }
     int numVars = vars.rows();
     VectorXd grad, x, H;
     DirectMetric s(schd.sDiagShift);
     double E0, stddev, rt;
     while (iter < maxIter)
     {
       E0 = 0.0;
       stddev = 0.0;
       rt = 0.0;
       grad.setZero(numVars);
       x.setZero(numVars + 1);
       H.setZero(numVars + 1);

       getMetric(vars, grad, H, s, E0, stddev, rt);
       write(vars);

       /*
         FOR DEBUGGING PURPOSE
       */
       /*
       s.BuildMetric();
       for(int i=1; i<vars.rows(); i++)
         s.S(i,i) += schd.sDiagShift;
       //MatrixXd s_inv = MatrixXd::Zero(vars.rows() + 1, vars.rows() + 1);
       //PInv(s.S,s_inv);
       //x = s_inv * H;
       */

       //xguess << 1.0, vars;
       //ConjGrad(s, s.H, xguess, schd.cgIter, x);

       x[0] = 1.0;
       ConjGrad(s, H, schd.cgIter, x);
       
       for (int i = 0; i < vars.rows(); i++)
       {
         vars(i) += (x(i+1) / x(0));
       }
       
#ifndef SERIAL
       MPI_Bcast(&(vars[0]), vars.rows(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
       if (commrank == 0)
         std::cout << format("%5i %14.8f (%8.2e) %14.8f %8.1f %10i %8.2f\n") % iter % E0 % stddev % (grad.norm()) % (rt) % (schd.stochasticIter) % ((getTime() - startofCalc));
       iter++;
     }
   }
};   

#endif
