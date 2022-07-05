//hello
#ifndef WalkerHelperMultiSlater_HEADER_H
#define WalkerHelperMultiSlater_HEADER_H

#include MultiSlater.h
#include WalkerHelperBase.h

template<>
class WalkerHelper<MultiSlater>
{

 public:
  MatrixXcd t;                    // A^{-1}
  complex<double> refOverlap;     // < n | phi_0 >
  std::vector<double> ciOverlaps; // Re (< n | phi_i >), include parity
  double totalOverlap;            // Re (< n | psi >)
  std::vector<complex<double>> ciOverlapRatios; // < n | phi_i > / < n | phi_0 >, include parity, for orb gradient
  complex<double> totalComplexOverlap;  // < n | psi >, for orb gradient
  std::vector<int> closedOrbs, openOrbs;    // set of closed and open orbitals in the walker
  MatrixXcd r, c, b, rt, tc, rtc_b;  // intermediate tables
  vector<Matrix4cd> tcSlice;      // intermediate tables

  WalkerHelper() {};
  
  WalkerHelper(const MultiSlater &w, const Determinant &d)
  {
    //fill the closed orbs for the walker
    closedOrbs.clear();
    openOrbs.clear();
    vector<int> closedBeta, openBeta;
    d.getOpenClosedAlphaBeta(openOrbs, closedOrbs, openBeta, closedBeta);
    for (int& c_i : closedBeta) c_i += Determinant::norbs;
    for (int& o_i : openBeta) o_i += Determinant::norbs;
    closedOrbs.insert(closedOrbs.end(), closedBeta.begin(), closedBeta.end());
    openOrbs.insert(openOrbs.end(), openBeta.begin(), openBeta.end());
    
    initInvDetsTables(w);
    if (commrank == 0) {
      if (refOverlap.real() == 0 || totalOverlap == 0) {
        cout << "refOverlap = " << refOverlap << ", totalOverlap = " << totalOverlap << endl;
        cout << "walker det: " << d << endl;
        cout << "ciOverlaps\n";
        for (double ci: ciOverlaps) cout << ci << endl;
      }
    }
  }

  void initInvDetsTables(const MultiSlater &w)
  {
    int norbs = Determinant::norbs;
  
    // inverse and refDet
    MatrixXcd a;
    Eigen::Map<VectorXi> occRows(&closedOrbs[0], closedOrbs.size());
    auto refCopy = w.ref;
    Eigen::Map<VectorXi> occColumns(&refCopy[0], refCopy.size());
    igl::slice(w.getHforbs(), occRows, occColumns, a); 
    Eigen::FullPivLU<MatrixXcd> lua(a);
    if (lua.isInvertible()) {
      t = lua.inverse();
      refOverlap = lua.determinant();
    }
    else {
      if (commrank == 0) {
        cout << "overlap with zeroth determinant not invertible" << endl;
        cout << "occRows\n" << occRows.transpose() <<  endl;
        cout << "occColumns\n" << occColumns.transpose() <<  endl;
        cout << "a\n" << a << endl;
        cout << "w.Hforbs\n" <<  w.Hforbs << endl;
      }
#ifndef SERIAL
  MPI_Barrier(MPI_COMM_WORLD);
#endif
      exit(0);
    }
    
    // tables
    // TODO: change table structure so that only unoccupied orbitals are present in r and c, 
    // this is not done currently because table updates are easier with all orbitals, but leads to bigger tables
    //VectorXi all = VectorXi::LinSpaced(2*norbs, 0, 2*norbs - 1);
    auto openCopy = w.open;
    Eigen::Map<VectorXi> openColumns(&openCopy[0], openCopy.size());
    Eigen::Map<VectorXi> openRows(&openOrbs[0], openOrbs.size());
    //igl::slice(w.getHforbs(), all, occColumns, r);
    //igl::slice(w.getHforbs(), occRows, all, c);
    igl::slice(w.getHforbs(), openRows, occColumns, r);
    igl::slice(w.getHforbs(), occRows, openColumns, c);
    igl::slice(w.getHforbs(), openRows, openColumns, b);
    rt = r * t;
    tc = t * c;
    rtc_b = rt * c - b;
    //rtc_b = rt * c - w.getHforbs();

    // overlaps with phi_i
    ciOverlaps.clear();
    ciOverlapRatios.clear();
    ciOverlaps.push_back(refOverlap.real());
    ciOverlapRatios.push_back(1.);
    totalOverlap = w.ciCoeffs[0] * ciOverlaps[0];
    totalComplexOverlap = w.ciCoeffs[0] * refOverlap;
    for (int i = 1; i < w.numDets; i++) {
      if (w.ciExcitations[i][0].size() == 4) {
        Matrix4cd sliceMat;
        igl::slice(tc, w.ciExcitations[i][0], w.ciExcitations[i][1], sliceMat);
        tcSlice.push_back(sliceMat);
        complex<double> ratio = sliceMat.determinant() * complex<double>(w.ciParity[i]);
        ciOverlapRatios.push_back(ratio);
        ciOverlaps.push_back((ratio * refOverlap).real());
        totalOverlap += w.ciCoeffs[i] * ciOverlaps[i];
        totalComplexOverlap += w.ciCoeffs[i] * ratio * refOverlap;
      }
      else {
        MatrixXcd sliceMat;
        igl::slice(tc, w.ciExcitations[i][0], w.ciExcitations[i][1], sliceMat);
        complex<double> ratio = calcDet(sliceMat) * complex<double>(w.ciParity[i]);
        ciOverlapRatios.push_back(ratio);
        ciOverlaps.push_back((ratio * refOverlap).real());
        totalOverlap += w.ciCoeffs[i] * ciOverlaps[i];
        totalComplexOverlap += w.ciCoeffs[i] * ratio * refOverlap;
      }
    }
  }

  void excitationUpdate(const MultiSlater &w, vector<int>& cre, vector<int>& des, bool sz, double parity, const Determinant& excitedDet)
  {
    // sherman morrison to update inverse
    // right now only table rt is updated efficiently
    auto refCopy = w.ref;
    Eigen::Map<VectorXi> occColumns(&refCopy[0], refCopy.size());
    MatrixXcd tOld = t;
    complex<double> overlapOld = refOverlap;
    MatrixXcd rtOld = rt;
    calculateInverseDeterminantWithRowChange(tOld, overlapOld, rtOld, t, refOverlap, rt, cre, des, occColumns, closedOrbs, w.getHforbs(), 0);
    refOverlap *= parity;
   
    int norbs = Determinant::norbs;
    closedOrbs.clear();
    openOrbs.clear();
    vector<int> closedBeta, openBeta;
    excitedDet.getOpenClosedAlphaBeta(openOrbs, closedOrbs, openBeta, closedBeta);
    for (int& c_i : closedBeta) c_i += norbs;
    for (int& o_i : openBeta) o_i += norbs;
    closedOrbs.insert(closedOrbs.end(), closedBeta.begin(), closedBeta.end());
    openOrbs.insert(openOrbs.end(), openBeta.begin(), openBeta.end());
    
    // TODO: these tables also need efficient updates
    auto openCopy = w.open;
    Eigen::Map<VectorXi> openColumns(&openCopy[0], openCopy.size());
    Eigen::Map<VectorXi> occRows(&closedOrbs[0], closedOrbs.size());
    Eigen::Map<VectorXi> openRows(&openOrbs[0], openOrbs.size());
    //VectorXi all = VectorXi::LinSpaced(2*norbs, 0, 2*norbs - 1);
    igl::slice(w.getHforbs(), openRows, occColumns, r);
    igl::slice(w.getHforbs(), occRows, openColumns, c);
    igl::slice(w.getHforbs(), openRows, openColumns, b);
    rt = r * t;
    tc = t * c;
    rtc_b = rt * c - b;
    //rtc_b = rt * c - w.getHforbs();
    
    // overlaps with phi_i
    ciOverlaps.clear();
    ciOverlapRatios.clear();
    ciOverlaps.push_back(refOverlap.real());
    ciOverlapRatios.push_back(1.);
    totalOverlap = w.ciCoeffs[0] * ciOverlaps[0];
    totalComplexOverlap = w.ciCoeffs[0] * refOverlap;
    size_t count4 = 0;
    for (int j = 1; j < w.numDets; j++) {
      int rank = w.ciExcitations[j][0].size();
      complex<double> sliceDet;
      if (rank == 1) sliceDet = tc(w.ciExcitations[j][0][0], w.ciExcitations[j][1][0]);
      
      else if (rank == 2) sliceDet = tc(w.ciExcitations[j][0][0], w.ciExcitations[j][1][0]) * tc(w.ciExcitations[j][0][1], w.ciExcitations[j][1][1])
- tc(w.ciExcitations[j][0][0], w.ciExcitations[j][1][1]) * tc(w.ciExcitations[j][0][1], w.ciExcitations[j][1][0]);
      
      else if (rank == 3) {
        //igl::slice(tc, w.ciExcitations[j][0], w.ciExcitations[j][1], sliceMat);
        Matrix3cd tcSlice;
        for  (int mu = 0; mu < 3; mu++) 
          for (int nu = 0; nu < 3; nu++) 
            tcSlice(mu, nu) = tc(w.ciExcitations[j][0][mu], w.ciExcitations[j][1][nu]);
        //sliceDet = calcDet(sliceMat);
        sliceDet = tcSlice.determinant();
      }
      
      else if (rank == 4) {
        //igl::slice(tc, w.ciExcitations[j][0], w.ciExcitations[j][1], sliceMat);
        for  (int mu = 0; mu < 4; mu++) 
          for (int nu = 0; nu < 4; nu++) 
            tcSlice[count4](mu, nu) = tc(w.ciExcitations[j][0][mu], w.ciExcitations[j][1][nu]);
        //sliceDet = calcDet(sliceMat);
        sliceDet = tcSlice[count4].determinant();
        count4++;
      }
      
      else {
        MatrixXcd sliceMat = MatrixXcd::Zero(rank, rank);
        //igl::slice(tc, w.ciExcitations[j][0], w.ciExcitations[j][1], sliceMat);
        for (int mu = 0; mu < rank; mu++) 
          for (int nu = 0; nu < rank; nu++) 
            sliceMat(mu, nu) = tc(w.ciExcitations[j][0][mu], w.ciExcitations[j][1][nu]);
        sliceDet = sliceMat.determinant();
      }
      
      ciOverlaps.push_back((sliceDet * refOverlap).real() * w.ciParity[j]);
      ciOverlapRatios.push_back(sliceDet * complex<double>(w.ciParity[j]));
      totalOverlap += w.ciCoeffs[j] * ciOverlaps[j];
      totalComplexOverlap += w.ciCoeffs[j] * ciOverlapRatios[j] * refOverlap;
    }
    
    if (commrank == 0) {
      if (refOverlap.real() == 0 || totalOverlap == 0) {
        cout << "refOverlap = " << refOverlap << ", totalOverlap = " << totalOverlap << endl;
        cout << "walker det: " << excitedDet << endl;
        cout << "ciOverlaps\n";
        for (double ci: ciOverlaps) cout << ci << endl;
      }
    }
  }

  void getRelIndices(int i, int &relI, int a, int &relA, bool sz) const
  {
    //relI = std::lower_bound(closedOrbs[sz].begin(), closedOrbs[sz].end(), i) - closedOrbs[sz].begin();
    relI = std::search_n(closedOrbs.begin(), closedOrbs.end(), 1, i + sz * Determinant::norbs) - closedOrbs.begin();
    relA = std::search_n(openOrbs.begin(), openOrbs.end(), 1, a + sz * Determinant::norbs) - openOrbs.begin();
    //relA = a + sz * Determinant::norbs;
  }

};

#endif
