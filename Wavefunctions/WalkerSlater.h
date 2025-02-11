#ifndef WalkerSlater_HEADER_H
#define WalkerSlater_HEADER_H
#include "Slater.h"
#include "WalkerBase.h"
#include "WalkerHelperSlater.h"
#include <unordered_set>
#include <boost/algorithm/string/predicate.hpp>


template<typename Corr>
struct Walker<Corr, Slater> {

  Determinant d;
  WalkerHelper<Corr> corrHelper;
  WalkerHelper<Slater> refHelper;
  unordered_set<int> excitedHoles;    //spin orbital indices of excited electrons (in core orbitals) in d
  unordered_set<int> excitedOrbs;     //spin orbital indices of excited electrons (in virtual orbitals) in d

  Walker() {};
  
  Walker(Corr &corr, const Slater &ref) 
  {
    initDet(ref.getHforbsA().real(), ref.getHforbsB().real());
    refHelper = WalkerHelper<Slater>(ref, d);
    corrHelper = WalkerHelper<Corr>(corr, d);
  }

  template<typename Wave>
  Walker(Wave &wave, const Determinant &pd)
  {
    d = pd;
    refHelper = WalkerHelper<Slater>(wave.ref, pd);
    corrHelper = WalkerHelper<Corr>(wave.corr, pd);
  }

  Walker(Corr &corr, const Slater &ref, const Determinant &pd) : d(pd), refHelper(ref, pd), corrHelper(corr, pd) {}; 

  Determinant& getDet() {return d;}
  void readBestDeterminant(Determinant& d) const 
  {
    if (commrank == 0) {
      char file[5000];
      sprintf(file, "BestDeterminant.txt");
      std::ifstream ifs(file, std::ios::binary);
      boost::archive::binary_iarchive load(ifs);
      load >> d;
    }
#ifndef SERIAL
    MPI_Bcast(&d.reprA, DetLen, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&d.reprB, DetLen, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
  }

  /**
   * makes det based on mo coeffs 
   */
  void guessBestDeterminant(Determinant& d, const Eigen::MatrixXd& HforbsA, const Eigen::MatrixXd& HforbsB) const 
  {
    int norbs = Determinant::norbs;
    int nalpha = Determinant::nalpha;
    int nbeta = Determinant::nbeta;

    d = Determinant();
    if (boost::iequals(schd.determinantFile, "")) {
      for (int i = 0; i < nalpha; i++) {
        int bestorb = 0;
        double maxovlp = 0;
        for (int j = 0; j < norbs; j++) {
          if (abs(HforbsA(i, j)) > maxovlp && !d.getoccA(j)) {
            maxovlp = abs(HforbsA(i, j));
            bestorb = j;
          }
        }
        d.setoccA(bestorb, true);
      }
      for (int i = 0; i < nbeta; i++) {
        int bestorb = 0;
        double maxovlp = 0;
        for (int j = 0; j < norbs; j++) {
          if (schd.hf == "rhf" || schd.hf == "uhf") {
            if (abs(HforbsB(i, j)) > maxovlp && !d.getoccB(j)) {
              bestorb = j;
              maxovlp = abs(HforbsB(i, j));
            }
          }
          else {
            if (abs(HforbsB(i+norbs, j)) > maxovlp && !d.getoccB(j)) {
              bestorb = j;
              maxovlp = abs(HforbsB(i+norbs, j));
            }
          }
        }
        d.setoccB(bestorb, true);
      }
    }
    else if (boost::iequals(schd.determinantFile, "bestDet")) {
      std::vector<Determinant> dets;
      std::vector<double> ci;
      readDeterminants(schd.determinantFile, dets, ci);
      d = dets[0];
    }
  }

  void initDet(const Eigen::MatrixXd& HforbsA, const Eigen::MatrixXd& HforbsB) 
  {
    bool readDeterminant = false;
    char file[5000];
    sprintf(file, "BestDeterminant.txt");

    {
      ifstream ofile(file);
      if (ofile)
        readDeterminant = true;
    }
    if (readDeterminant)
      readBestDeterminant(d);
    else
      guessBestDeterminant(d, HforbsA, HforbsB);
  }

  double getIndividualDetOverlap(int i) const
  {
    return (refHelper.thetaDet[i][0] * refHelper.thetaDet[i][1]).real();
  }

  double getDetOverlap(const Slater &ref) const
  {
    double ovlp = 0.0;
    for (int i = 0; i < refHelper.thetaDet.size(); i++) {
      ovlp += ref.getciExpansion()[i] * (refHelper.thetaDet[i][0] * refHelper.thetaDet[i][1]).real();
    }
    return ovlp;
  }

  double getDetFactor(int i, int a, const Slater &ref) const 
  {
    if (i % 2 == 0)
      return getDetFactor(i / 2, a / 2, 0, ref);
    else                                   
      return getDetFactor(i / 2, a / 2, 1, ref);
  }

  double getDetFactor(int I, int J, int A, int B, const Slater &ref) const 
  {
    if (I % 2 == J % 2 && I % 2 == 0)
      return getDetFactor(I / 2, J / 2, A / 2, B / 2, 0, 0, ref);
    else if (I % 2 == J % 2 && I % 2 == 1)                  
      return getDetFactor(I / 2, J / 2, A / 2, B / 2, 1, 1, ref);
    else if (I % 2 != J % 2 && I % 2 == 0)                  
      return getDetFactor(I / 2, J / 2, A / 2, B / 2, 0, 1, ref);
    else                                                    
      return getDetFactor(I / 2, J / 2, A / 2, B / 2, 1, 0, ref);
  }

  double getDetFactor(int i, int a, bool sz, const Slater &ref) const
  {
    int tableIndexi, tableIndexa;
    refHelper.getRelIndices(i, tableIndexi, a, tableIndexa, sz); 

    double detFactorNum = 0.0;
    double detFactorDen = 0.0;
    for (int j = 0; j < ref.getDeterminants().size(); j++)
    {
      double factor = (refHelper.rTable[j][sz](tableIndexa, tableIndexi) * refHelper.thetaDet[j][0] * refHelper.thetaDet[j][1]).real() * ref.getciExpansion()[j] /  getDetOverlap(ref);
      detFactorNum += ref.getciExpansion()[j] * factor * (refHelper.thetaDet[j][0] * refHelper.thetaDet[j][1]).real();
      detFactorDen += ref.getciExpansion()[j] * (refHelper.thetaDet[j][0] * refHelper.thetaDet[j][1]).real();
    }
    return detFactorNum / detFactorDen;
  }
  
  double getDetFactor(int i, int j, int a, int b, bool sz1, bool sz2, const Slater &ref) const
  {
    int tableIndexi, tableIndexa, tableIndexj, tableIndexb;
    refHelper.getRelIndices(i, tableIndexi, a, tableIndexa, sz1); 
    refHelper.getRelIndices(j, tableIndexj, b, tableIndexb, sz2) ;

    double detFactorNum = 0.0;
    double detFactorDen = 0.0;
    for (int j = 0; j < ref.getDeterminants().size(); j++)
    {
      double factor;
      if (sz1 == sz2 || refHelper.hftype == 2)
        factor =((refHelper.rTable[j][sz1](tableIndexa, tableIndexi) * refHelper.rTable[j][sz1](tableIndexb, tableIndexj) 
            - refHelper.rTable[j][sz1](tableIndexb, tableIndexi) *refHelper.rTable[j][sz1](tableIndexa, tableIndexj)) * refHelper.thetaDet[j][0] * refHelper.thetaDet[j][1]).real() * ref.getciExpansion()[j]/ getDetOverlap(ref);
      else
        factor = (refHelper.rTable[j][sz1](tableIndexa, tableIndexi) * refHelper.rTable[j][sz2](tableIndexb, tableIndexj) * refHelper.thetaDet[j][0] * refHelper.thetaDet[j][1]).real() * ref.getciExpansion()[j]/ getDetOverlap(ref);
      detFactorNum += ref.getciExpansion()[j] * factor * (refHelper.thetaDet[j][0] * refHelper.thetaDet[j][1]).real();
      detFactorDen += ref.getciExpansion()[j] * (refHelper.thetaDet[j][0] * refHelper.thetaDet[j][1]).real();
    }
    return detFactorNum / detFactorDen;
  }
 
  //only works for ghf
  double getDetFactor(std::array<unordered_set<int>, 2> &from, std::array<unordered_set<int>, 2> &to) const
  {
    if (from[0].size() + from[1].size() == 0) return 1.;
    int numExc = from[0].size() + from[1].size();
    Eigen::VectorXi tableIndicesRow = Eigen::VectorXi::Zero(from[0].size() + from[1].size());
    Eigen::VectorXi tableIndicesCol = Eigen::VectorXi::Zero(from[0].size() + from[1].size());
    Determinant dcopy = d;
    double parity = 1.;
    int count = 0;
    for (int sz = 0; sz < 2; sz++) {//iterate over spins
      auto itFrom = from[sz].begin();
      auto itTo = to[sz].begin();
      for (int n = 0; n < from[sz].size(); n++) {//iterate over excitations
        int i = *itFrom, a = *itTo;
        itFrom = std::next(itFrom); itTo = std::next(itTo);
        refHelper.getRelIndices(i, tableIndicesCol(count), a, tableIndicesRow(count), sz);
        count++;
        parity *= dcopy.parity(a, i, sz);
        dcopy.setocc(i, sz, false);  
        dcopy.setocc(a, sz, true);
      }
    }

    Eigen::MatrixXcd detSlice = Eigen::MatrixXcd::Zero(numExc,numExc);
    igl::slice(refHelper.rTable[0][0], tableIndicesRow, tableIndicesCol, detSlice);
    complex<double> det(0.,0.);
    if (detSlice.rows() == 1) det = detSlice(0, 0);
    else if (detSlice.rows() == 2) det = detSlice(0, 0) * detSlice(1, 1) - detSlice(0, 1) * detSlice(1, 0);
    else if (detSlice.rows() == 3) det =   detSlice(0, 0) * (detSlice(1, 1) * detSlice(2, 2) - detSlice(1, 2) * detSlice(2, 1))
                                         - detSlice(0, 1) * (detSlice(1, 0) * detSlice(2, 2) - detSlice(1, 2) * detSlice(2, 0))
                                         + detSlice(0, 2) * (detSlice(1, 0) * detSlice(2, 1) - detSlice(1, 1) * detSlice(2, 0));
    else if (detSlice.rows() == 4) det = detSlice(0,3) * detSlice(1,2) * detSlice(2,1) * detSlice(3,0) - detSlice(0,2) * detSlice(1,3) * detSlice(2,1) * detSlice(3,0) -
       detSlice(0,3) * detSlice(1,1) * detSlice(2,2) * detSlice(3,0) + detSlice(0,1) * detSlice(1,3) * detSlice(2,2) * detSlice(3,0) +
       detSlice(0,2) * detSlice(1,1) * detSlice(2,3) * detSlice(3,0) - detSlice(0,1) * detSlice(1,2) * detSlice(2,3) * detSlice(3,0) -
       detSlice(0,3) * detSlice(1,2) * detSlice(2,0) * detSlice(3,1) + detSlice(0,2) * detSlice(1,3) * detSlice(2,0) * detSlice(3,1) +
       detSlice(0,3) * detSlice(1,0) * detSlice(2,2) * detSlice(3,1) - detSlice(0,0) * detSlice(1,3) * detSlice(2,2) * detSlice(3,1) -
       detSlice(0,2) * detSlice(1,0) * detSlice(2,3) * detSlice(3,1) + detSlice(0,0) * detSlice(1,2) * detSlice(2,3) * detSlice(3,1) +
       detSlice(0,3) * detSlice(1,1) * detSlice(2,0) * detSlice(3,2) - detSlice(0,1) * detSlice(1,3) * detSlice(2,0) * detSlice(3,2) -
       detSlice(0,3) * detSlice(1,0) * detSlice(2,1) * detSlice(3,2) + detSlice(0,0) * detSlice(1,3) * detSlice(2,1) * detSlice(3,2) +
       detSlice(0,1) * detSlice(1,0) * detSlice(2,3) * detSlice(3,2) - detSlice(0,0) * detSlice(1,1) * detSlice(2,3) * detSlice(3,2) -
       detSlice(0,2) * detSlice(1,1) * detSlice(2,0) * detSlice(3,3) + detSlice(0,1) * detSlice(1,2) * detSlice(2,0) * detSlice(3,3) +
       detSlice(0,2) * detSlice(1,0) * detSlice(2,1) * detSlice(3,3) - detSlice(0,0) * detSlice(1,2) * detSlice(2,1) * detSlice(3,3) -
       detSlice(0,1) * detSlice(1,0) * detSlice(2,2) * detSlice(3,3) + detSlice(0,0) * detSlice(1,1) * detSlice(2,2) * detSlice(3,3);

    //complex<double> det = detSlice.determinant();
    double num = (det * refHelper.thetaDet[0][0] * refHelper.thetaDet[0][1]).real();
    double den = (refHelper.thetaDet[0][0] * refHelper.thetaDet[0][1]).real();
    return parity * num / den;
  }

  void update(int i, int a, bool sz, const Slater &ref, Corr &corr, bool doparity = true)
  {
    double p = 1.0;
    if (doparity) p *= d.parity(a, i, sz);
    d.setocc(i, sz, false);
    d.setocc(a, sz, true);
    if (refHelper.hftype == Generalized) {
      int norbs = Determinant::norbs;
      vector<int> cre{ a + sz * norbs }, des{ i + sz * norbs };
      refHelper.excitationUpdateGhf(ref, cre, des, sz, p, d);
    }
    else
    {
      vector<int> cre{ a }, des{ i };
      refHelper.excitationUpdate(ref, cre, des, sz, p, d);
    }

    corrHelper.updateHelper(corr, d, i, a, sz);
  }

  void update(int i, int j, int a, int b, bool sz, const Slater &ref, Corr& corr, bool doparity = true)
  {
    double p = 1.0;
    Determinant dcopy = d;
    if (doparity) p *= d.parity(a, i, sz);
    d.setocc(i, sz, false);
    d.setocc(a, sz, true);
    if (doparity) p *= d.parity(b, j, sz);
    d.setocc(j, sz, false);
    d.setocc(b, sz, true);
    if (refHelper.hftype == Generalized) {
      int norbs = Determinant::norbs;
      vector<int> cre{ a + sz * norbs, b + sz * norbs }, des{ i + sz * norbs, j + sz * norbs };
      refHelper.excitationUpdateGhf(ref, cre, des, sz, p, d);
    }
    else {
      vector<int> cre{ a, b }, des{ i, j };
      refHelper.excitationUpdate(ref, cre, des, sz, p, d);
    }
    corrHelper.updateHelper(corr, d, i, j, a, b, sz);
  }

  void updateWalker(const Slater &ref, Corr& corr, int ex1, int ex2, bool doparity = true)
  {
    int norbs = Determinant::norbs;
    int I = ex1 / 2 / norbs, A = ex1 - 2 * norbs * I;
    int J = ex2 / 2 / norbs, B = ex2 - 2 * norbs * J;
    if (I % 2 == J % 2 && ex2 != 0) {
      if (I % 2 == 1) {
        update(I / 2, J / 2, A / 2, B / 2, 1, ref, corr, doparity);
      }
      else {
        update(I / 2, J / 2, A / 2, B / 2, 0, ref, corr, doparity);
      }
    }
    else {
      if (I % 2 == 0)
        update(I / 2, A / 2, 0, ref, corr, doparity);
      else
        update(I / 2, A / 2, 1, ref, corr, doparity);

      if (ex2 != 0) {
        if (J % 2 == 1) {
          update(J / 2, B / 2, 1, ref, corr, doparity);
        }
        else {
          update(J / 2, B / 2, 0, ref, corr, doparity);
        }
      }
    }
  }

  void exciteWalker(const Slater &ref, Corr& corr, int excite1, int excite2, int norbs)
  {
    int I1 = excite1 / (2 * norbs), A1 = excite1 % (2 * norbs);

    if (I1 % 2 == 0)
      update(I1 / 2, A1 / 2, 0, ref, corr);
    else
      update(I1 / 2, A1 / 2, 1, ref, corr);

    if (excite2 != 0) {
      int I2 = excite2 / (2 * norbs), A2 = excite2 % (2 * norbs);
      if (I2 % 2 == 0)
        update(I2 / 2, A2 / 2, 0, ref, corr);
      else
        update(I2 / 2, A2 / 2, 1, ref, corr);
    }
  }

  void OverlapWithOrbGradient(const Slater &ref, Eigen::VectorXd &grad, double detovlp) const
  {
    int norbs = Determinant::norbs;
    Determinant walkerDet = d;

    //K and L are relative row and col indices
    int KA = 0, KB = 0;
    for (int k = 0; k < norbs; k++) { //walker indices on the row
      if (walkerDet.getoccA(k)) {
        for (int det = 0; det < ref.getDeterminants().size(); det++) {
          Determinant refDet = ref.getDeterminants()[det];
          int L = 0;
          for (int l = 0; l < norbs; l++) {
            if (refDet.getoccA(l)) {
              grad(2 * k * norbs + 2 * l) += ref.getciExpansion()[det] * (refHelper.thetaInv[0](L, KA) * refHelper.thetaDet[det][0] * refHelper.thetaDet[det][1]).real() /detovlp;
              grad(2 * k * norbs + 2 * l + 1) += ref.getciExpansion()[det] * (- refHelper.thetaInv[0](L, KA) * refHelper.thetaDet[det][0] * refHelper.thetaDet[det][1]).imag() /detovlp;
              L++;
            }
          }
        }
        KA++;
      }
      if (walkerDet.getoccB(k)) {
        for (int det = 0; det < ref.getDeterminants().size(); det++) {
          Determinant refDet = ref.getDeterminants()[det];
          int L = 0;
          for (int l = 0; l < norbs; l++) {
            if (refDet.getoccB(l)) {
              if (refHelper.hftype == UnRestricted) {
                grad(2 * norbs * norbs + 2 * k * norbs + 2 * l) += ref.getciExpansion()[det] * (refHelper.thetaInv[1](L, KB) * refHelper.thetaDet[det][0] * refHelper.thetaDet[det][1]).real() / detovlp;
                grad(2 * norbs * norbs + 2 * k * norbs + 2 * l + 1) += ref.getciExpansion()[det] * (- refHelper.thetaInv[1](L, KB) * refHelper.thetaDet[det][0] * refHelper.thetaDet[det][1]).imag() / detovlp;
              }
              else {
                grad(2 * k * norbs + 2 * l) += ref.getciExpansion()[det] * (refHelper.thetaInv[1](L, KB) * refHelper.thetaDet[det][0] * refHelper.thetaDet[det][1]).real() / detovlp;
                grad(2 * k * norbs + 2 * l + 1) += ref.getciExpansion()[det] * (- refHelper.thetaInv[1](L, KB) * refHelper.thetaDet[det][0] * refHelper.thetaDet[det][1]).imag() / detovlp;
              }
              L++;
            }
          }
        }
        KB++;
      }
    }
  }

  void OverlapWithOrbGradientGhf(const Slater &ref, Eigen::VectorXd &grad, double detovlp) const
  {
    int norbs = Determinant::norbs;
    Determinant walkerDet = d;
    Determinant refDet = ref.getDeterminants()[0];

    //K and L are relative row and col indices
    int K = 0;
    for (int k = 0; k < norbs; k++) { //walker indices on the row
      if (walkerDet.getoccA(k)) {
        int L = 0;
        for (int l = 0; l < norbs; l++) {
          if (refDet.getoccA(l)) {
            grad(4 * k * norbs + 2 * l) += (refHelper.thetaInv[0](L, K) * refHelper.thetaDet[0][0]).real() / detovlp;
            grad(4 * k * norbs + 2 * l + 1) += (- refHelper.thetaInv[0](L, K) * refHelper.thetaDet[0][0]).imag() / detovlp;
            L++;
          }
        }
        for (int l = 0; l < norbs; l++) {
          if (refDet.getoccB(l)) {
            grad(4 * k * norbs + 2 * norbs + 2 * l) += (refHelper.thetaInv[0](L, K) * refHelper.thetaDet[0][0]).real() / detovlp;
            grad(4 * k * norbs + 2 * norbs + 2 * l + 1) += (- refHelper.thetaInv[0](L, K) * refHelper.thetaDet[0][0]).imag() / detovlp;
            L++;
          }
        }
        K++;
      }
    }
    for (int k = 0; k < norbs; k++) { //walker indices on the row
      if (walkerDet.getoccB(k)) {
        int L = 0;
        for (int l = 0; l < norbs; l++) {
          if (refDet.getoccA(l)) {
            grad(4 * norbs * norbs +  4 * k * norbs + 2 * l) += (refHelper.thetaDet[0][0] * refHelper.thetaInv[0](L, K)).real() / detovlp;
            grad(4 * norbs * norbs +  4 * k * norbs + 2 * l + 1) += (- refHelper.thetaDet[0][0] * refHelper.thetaInv[0](L, K)).imag() / detovlp;
            L++;
          } 
        }
        for (int l = 0; l < norbs; l++) {
          if (refDet.getoccB(l)) {
            grad(4 * norbs * norbs +  4 * k * norbs + 2 * norbs + 2 * l) += (refHelper.thetaDet[0][0] * refHelper.thetaInv[0](L, K)).real() / detovlp;
            grad(4 * norbs * norbs +  4 * k * norbs + 2 * norbs + 2 * l + 1) += (- refHelper.thetaDet[0][0] * refHelper.thetaInv[0](L, K)).imag() / detovlp;
            L++;
          }
        }
        K++;
      }
    }
  }

  void OverlapWithGradient(const Slater &ref, Eigen::VectorBlock<Eigen::VectorXd> &grad) const
  {
    double detovlp = getDetOverlap(ref);
    //for (int i = 0; i < ref.ciExpansion.size(); i++)
    //  grad[i] += getIndividualDetOverlap(i) / detovlp;
    grad[0] = 0.;
    if (ref.determinants.size() <= 1 && schd.optimizeOrbs) {
      //if (hftype == UnRestricted)
      Eigen::VectorXd gradOrbitals;
      if (ref.hftype == UnRestricted) {
        gradOrbitals = Eigen::VectorXd::Zero(4 * ref.HforbsA.rows() * ref.HforbsA.rows());
        OverlapWithOrbGradient(ref, gradOrbitals, detovlp);
      }
      else {
        gradOrbitals = Eigen::VectorXd::Zero(2 * ref.HforbsA.rows() * ref.HforbsA.rows());
        if (ref.hftype == Restricted) OverlapWithOrbGradient(ref, gradOrbitals, detovlp);
        else OverlapWithOrbGradientGhf(ref, gradOrbitals, detovlp);
      }
      for (int i = 0; i < gradOrbitals.size(); i++)
        grad[ref.ciExpansion.size() + i] += gradOrbitals[i];
    }
    //cout << "ref grad\n" << grad << endl;
  }

  friend ostream& operator<<(ostream& os, const Walker<Corr, Slater>& w) {
    os << w.d << endl << endl;
    os << "alphaTable\n" << w.refHelper.rTable[0][0] << endl << endl;
    os << "betaTable\n" << w.refHelper.rTable[0][1] << endl << endl;
    os << "dets\n" << w.refHelper.thetaDet[0][0] << "  " << w.refHelper.thetaDet[0][1] << endl << endl;
    os << "alphaInv\n" << w.refHelper.thetaInv[0] << endl << endl;
    os << "betaInv\n" << w.refHelper.thetaInv[1] << endl << endl;
    return os;
  }

};


#endif
