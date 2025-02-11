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

#ifndef SJastrow_HEADER_H
#define SJastrow_HEADER_H
#include "Correlator.h"
#include <Eigen/Dense>
#include <vector>
#include <iostream>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <string>

class Determinant;

/*
 * SJastrow is a product of correlators 
 */
class SJastrow {
 private:
  friend class boost::serialization::access;
  template<class Archive>
  void serialize (Archive & ar, const unsigned int version) {
    ar & SpinCorrelator;
  }
 public:
  Eigen::MatrixXd SpinCorrelator; //only lower triangle is relevant
  
  //reads correlator file and makes cpsArray, orbitalToSJastrow
  SJastrow ();
  
  double Overlap(const Determinant& d) const ;

  const double& operator()(int i, int j) const {return SpinCorrelator(std::max(i,j), std::min(i,j));}
  
  /*
   * Takes an occupation number representation of two determinants
   * in the local orbital basis and calculates the ratio of overlaps 
   * <d1|SJastrow>/<d2|SJastrow>
   * PARAMS:
   * 
   * Determinant: the occupation number representation of dets
   * 
   * RETURN:
   * <d1|SJastrow>/<d2|SJastrow>
   *
   */
  double OverlapRatio(const Determinant& d1, const Determinant& d2) const ;
  
  /*
   * return ratio of overlaps of SJastrow with d and (i->a,j->b)excited-d (=dcopy)
   */
  double OverlapRatio(int i, int a, const Determinant &dcopy, const Determinant &d) const ;
  double OverlapRatio(int i, int j, int a, int b, const Determinant &dcopy, const Determinant &d) const;
  
  /*
   * Takes an occupation number representation of a determinant
   * in the local orbital basis and calculates the overlap 
   * the SJastrow and also the overlap of the determinant
   * with respect to the tangent vector w.r.t to all the 
   * parameters in the wavefuncion
   * 
   * PARAMS:
   * 
   * Determinant: the occupation number representation of the determinant
   * grad       : the vector of the gradient. This vector is long and contains
   *              the space for storing gradients with respect to all the parameters
   *              in the wavefunction and not just this SJastrow and so the index startIndex
   *              is needed to indentify at what index should we start populating the
   *              gradient w.r.t to this SJastrow in the vector
   * ovlp       : <d|SJastrow>
   * startIndex : The location in the vector gradient from where to start populating
   *              the gradient w.r.t to the current Correlator
   */
  void   OverlapWithGradient  (const Determinant& d, 
                               Eigen::VectorBlock<Eigen::VectorXd>& grad,
                               const double& ovlp) const;

  void getVariables(Eigen::VectorBlock<Eigen::VectorXd> &v) const ;
  long getNumVariables() const;
  void updateVariables(const Eigen::VectorBlock<Eigen::VectorXd> &v);
  void printVariables() const;
  void addNoise() ;
  std::string getfileName() const {return "SJastrow";};
};


#endif
