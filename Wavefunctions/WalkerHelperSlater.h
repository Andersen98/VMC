#ifndef WalkerHelperSlater_HEADER_H
#define WalkerHelperSlater_HEADER_H
#include "Slater.h"
#include "WalkerHelperBase.h"
#include "ShermanMorrisonWoodbury.h"
#include <igl/slice.h>

template<>
class WalkerHelper<Slater>
{

 public:
  HartreeFock hftype;                           //hftype same as that in slater
  std::array<Eigen::MatrixXcd, 2> thetaInv;          //inverse of the theta matrix
  std::vector<std::array<complex<double>, 2>> thetaDet;    //determinant of the theta matrix, vector for multidet
  std::array<vector<int>, 2> openOrbs;       //set of open orbitals in the walker
  std::array<vector<int>, 2> closedOrbs;     //set of closed orbitals in the walker
  std::array<vector<int>, 2> closedOrbsRef;  //set of closed orbitals in the reference (zeroth det)
  std::vector<std::array<Eigen::MatrixXcd, 2>> rTable;    //table used for efficiently, vector for multidet

  WalkerHelper() {};
  
  WalkerHelper(const Slater &w, const Determinant &d) 
  {
    hftype = w.hftype;
 
    //fill the spin strings for the walker and the zeroth reference det
    fillOpenClosedOrbs(d);
    closedOrbsRef[0].clear();
    closedOrbsRef[1].clear();
    w.getDeterminants()[0].getClosedAlphaBeta(closedOrbsRef[0], closedOrbsRef[1]);
  
    rTable.resize(w.getNumOfDets());
    thetaDet.resize(w.getNumOfDets());
  
    if (hftype == Generalized) {
      initInvDetsTablesGhf(w);
    }
    else {
      initInvDetsTables(w);
    }
  }

  void fillOpenClosedOrbs(const Determinant &d)
  {
    openOrbs[0].clear();
    openOrbs[1].clear();
    closedOrbs[0].clear();
    closedOrbs[1].clear();
    d.getOpenClosedAlphaBeta(openOrbs[0], closedOrbs[0], openOrbs[1], closedOrbs[1]);
  }

  void makeTable(const Slater &w, const Eigen::MatrixXcd& inv, const Eigen::Map<Eigen::VectorXi>& colClosed, int detIndex, bool sz)
  {
    Eigen::Map<Eigen::VectorXi> rowOpen(&openOrbs[sz][0], openOrbs[sz].size());
    rTable[detIndex][sz] = Eigen::MatrixXcd::Zero(openOrbs[sz].size(), closedOrbs[sz].size()); 
    Eigen::MatrixXcd HfopenTheta;
    igl::slice(w.getHforbs(sz), rowOpen, colClosed, HfopenTheta);
    rTable[detIndex][sz] = HfopenTheta * inv;
  }

  //void calcOtherDetsTables(const Slater& w, bool sz)
  //{
  //  Eigen::Map<Eigen::VectorXi> rowClosed(&closedOrbs[sz][0], closedOrbs[sz].size());
  //  vector<int> cre(closedOrbs[sz].size(), -1), des(closedOrbs[sz].size(), -1);
  //  for (int x = 1; x < w.getNumOfDets(); x++) {
  //    Eigen::MatrixXcd invCurrent;
  //    vector<int> ref;
  //    w.getDeterminants()[x].getClosed(sz, ref);
  //    Eigen::Map<Eigen::VectorXi> colClosed(&ref[0], ref.size());
  //    getDifferenceInOccupation(w.getDeterminants()[x], w.getDeterminants()[0], cre, des, sz);
  //    double parity = w.getDeterminants()[0].parity(cre, des, sz);
  //    calculateInverseDeterminantWithColumnChange(thetaInv[sz], thetaDet[0][sz], invCurrent, thetaDet[x][sz], cre, des, rowClosed, closedOrbsRef[sz], w.getHforbs(sz));
  //    thetaDet[x][sz] *= parity;
  //    makeTable(w, invCurrent, colClosed, x, sz);
  //  }
  //}

  //commenting out calcotherdetstables, uncomment for multidet
  void initInvDetsTables(const Slater &w)
  {
    int norbs = Determinant::norbs;
    for (int sz = 0; sz < 2; sz++) {
      Eigen::Map<Eigen::VectorXi> rowClosed(&closedOrbs[sz][0], closedOrbs[sz].size());
      Eigen::Map<Eigen::VectorXi> colClosed(&closedOrbsRef[sz][0], closedOrbsRef[sz].size());
      Eigen::MatrixXcd theta;
      igl::slice(w.getHforbs(sz), rowClosed, colClosed, theta); 
      Eigen::FullPivLU<Eigen::MatrixXcd> lua(theta);
      if (lua.isInvertible()) {
        thetaInv[sz] = lua.inverse();
        thetaDet[0][sz] = lua.determinant();
      }
      else {
        cout << sz << " overlap with determinant not invertible" << endl;
        exit(0);
      }
      rTable[0][sz] = Eigen::MatrixXcd::Zero(norbs, closedOrbs[sz].size()); 
      rTable[0][sz] = w.getHforbs(sz).block(0, 0, norbs, closedOrbs[sz].size()) * thetaInv[sz];
      
      //makeTable(w, thetaInv[sz], colClosed, 0, sz);
      //calcOtherDetsTables(w, sz);
    }
  }

  void concatenateGhf(const vector<int>& v1, const vector<int>& v2, vector<int>& result) const
  {
    int norbs = Determinant::norbs;
    result.clear();
    result = v1;
    result.insert(result.end(), v2.begin(), v2.end());    
    for (int j = v1.size(); j < v1.size() + v2.size(); j++)
      result[j] += norbs;
  }

  void makeTableGhf(const Slater &w, const Eigen::Map<Eigen::VectorXi>& colTheta)
  {
    rTable[0][0] = Eigen::MatrixXcd::Zero(openOrbs[0].size() + openOrbs[1].size(), closedOrbs[0].size() + closedOrbs[1].size()); 
    Eigen::MatrixXcd ghfOpen;
    vector<int> rowVec;
    concatenateGhf(openOrbs[0], openOrbs[1], rowVec);
    Eigen::Map<Eigen::VectorXi> rowOpen(&rowVec[0], rowVec.size());
    igl::slice(w.getHforbs(), rowOpen, colTheta, ghfOpen);
    rTable[0][0] = ghfOpen * thetaInv[0];
    rTable[0][1] = rTable[0][0];
  }

  void initInvDetsTablesGhf(const Slater &w)
  {
    int norbs = Determinant::norbs;
    vector<int> workingVec0, workingVec1;
    concatenateGhf(closedOrbs[0], closedOrbs[1], workingVec0);
    Eigen::Map<Eigen::VectorXi> rowTheta(&workingVec0[0], workingVec0.size());
    concatenateGhf(closedOrbsRef[0], closedOrbsRef[1], workingVec1);
    Eigen::Map<Eigen::VectorXi> colTheta(&workingVec1[0], workingVec1.size());
  
    Eigen::MatrixXcd theta;
    igl::slice(w.getHforbs(), rowTheta, colTheta, theta); 
    Eigen::FullPivLU<Eigen::MatrixXcd> lua(theta);
    if (lua.isInvertible()) {
      thetaInv[0] = lua.inverse();
      thetaDet[0][0] = lua.determinant();
    }
    else {
      Eigen::Map<Eigen::VectorXi> v1(&closedOrbs[0][0], closedOrbs[0].size());
      Eigen::Map<Eigen::VectorXi> v2(&closedOrbs[1][0], closedOrbs[1].size());
      cout << "alphaClosed\n" << v1 << endl << endl;
      cout << "betaClosed\n" << v2 << endl << endl;
      cout << "col\n" << colTheta << endl << endl;
      cout << theta << endl << endl;
      cout << "overlap with theta determinant not invertible" << endl;
      exit(0);
    }
    thetaDet[0][1] = 1.;
    rTable[0][0] = Eigen::MatrixXcd::Zero(2*norbs, closedOrbs[0].size() + closedOrbs[1].size()); 
    rTable[0][0] = w.getHforbs().block(0, 0, 2*norbs, closedOrbs[0].size() + closedOrbs[1].size()) * thetaInv[0];
    rTable[0][1] = rTable[0][0];
    //makeTableGhf(w, colTheta);
  }

  //commenting out calcotherdetstables, uncomment for multidet
  void excitationUpdate(const Slater &w, vector<int>& cre, vector<int>& des, bool sz, double parity, const Determinant& excitedDet)
  {
    Eigen::MatrixXcd invOld = thetaInv[sz];
    std::complex<double> detOld = thetaDet[0][sz];
    Eigen::MatrixXcd tableOld = rTable[0][sz];
    Eigen::Map<Eigen::VectorXi> colClosed(&closedOrbsRef[sz][0], closedOrbsRef[sz].size());
    calculateInverseDeterminantWithRowChange(invOld, detOld, tableOld, thetaInv[sz], thetaDet[0][sz], rTable[0][sz], cre, des, colClosed, closedOrbs[sz], w.getHforbs(sz), 1);
    thetaDet[0][sz] *= parity;
    fillOpenClosedOrbs(excitedDet);
    //makeTable(w, thetaInv[sz], colClosed, 0, sz);
    //calcOtherDetsTables(w, sz);
  }

  void excitationUpdateGhf(const Slater &w, vector<int>& cre, vector<int>& des, bool sz, double parity, const Determinant& excitedDet)
  {
    vector<int> colVec;
    concatenateGhf(closedOrbsRef[0], closedOrbsRef[1], colVec);
    Eigen::Map<Eigen::VectorXi> colTheta(&colVec[0], colVec.size());
    vector<int> rowIn;
    concatenateGhf(closedOrbs[0], closedOrbs[1], rowIn);
    Eigen::MatrixXcd invOld = thetaInv[0];
    std::complex<double> detOld = thetaDet[0][0];
    Eigen::MatrixXcd tableOld = rTable[0][0];
    calculateInverseDeterminantWithRowChange(invOld, detOld, tableOld, thetaInv[0], thetaDet[0][0], rTable[0][0], cre, des, colTheta, rowIn, w.getHforbs(), 1);
    rTable[0][1] = rTable[0][0];
    thetaDet[0][0] *= parity;
    fillOpenClosedOrbs(excitedDet);
    //makeTableGhf(w, colTheta);
  }

  void getRelIndices(int i, int &relI, int a, int &relA, bool sz) const 
  {
    //relI = std::lower_bound(closedOrbs[sz].begin(), closedOrbs[sz].end(), i) - closedOrbs[sz].begin();
    //relA = std::lower_bound(openOrbs[sz].begin(), openOrbs[sz].end(), a) - openOrbs[sz].begin();
    int factor = 0;
    if (hftype == 2 && sz != 0) factor = 1;
    relI = std::search_n(closedOrbs[sz].begin(), closedOrbs[sz].end(), 1, i) - closedOrbs[sz].begin() + factor * closedOrbs[0].size();
    //relA = std::search_n(openOrbs[sz].begin(), openOrbs[sz].end(), 1, a) - openOrbs[sz].begin() + factor * openOrbs[0].size();
    relA = a + factor * Determinant::norbs;
  }

};

#endif
