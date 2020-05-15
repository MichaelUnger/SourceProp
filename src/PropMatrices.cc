#include "PropMatrices.h"

#include <iostream>
#include <sstream>
#include <stdexcept>

#include<TClass.h>
#include <TFile.h>
#include <TKey.h>
#include <TH2D.h>

#include <utl/Units.h>
using namespace std;

namespace prop {

  // if these vectors are modified they must be initialized so that they are
  // in ascending order numerically
  const std::vector<double> DminGrid = {0, 1, 5, 7, 10, 20, 30, 40, 50, 75, 100};
  const double DminMin = DminGrid.front();
  const double DminMax = DminGrid.back();
  const double dM = 0.2;
  const double dZ0 = 0.25;
  const double Mmin = -5.0;
  const double Mmax = 5.0;
  const double Z0min = 0.;
  const double Z0max = 5.0;
  
  PropMatrices::PropMatrices(const double lgEmin,
                             const double lgEmax):
    fLgEmin(lgEmin), fLgEmax(lgEmax), fMaxDistance(0)
  {}

  bool
  PropMatrices::HasPrimary(const int Aprim)
    const
  {
    return fMatrices.find(Aprim) != fMatrices.end();
  }

  bool
  PropMatrices::HasMatrix(const int Aprim,
                          const int Asec)
    const
  {
    auto m = fMatrices.find(Aprim);
    if (m != fMatrices.end())
      return (m->second.find(Asec) != m->second.end());
    else
      return false;
  }



  TMatrixD&
  PropMatrices::GetMatrix(const int Aprim,
                          const int Asec)
  {
    return fMatrices[Aprim][Asec];
  }

  void
  PropMatrices::ResetPrimaryMap()
  {
    fMatrices.clear();
    return;
  }

  unsigned int
  PropMatrices::GetN()
    const
  {
    // size of first matrix in map
    for (auto& iter1 : fMatrices) {
      for (auto& iter2 : iter1.second) {
        const TMatrixD& m = iter2.second;
        if (m.GetNcols() == m.GetNrows())
          return m.GetNcols();
        else {
          cerr << " PropMatrices::GetN() -- Error: nCol != nRow?? " << endl;
          return 0;
        }
      }
    }
    return 0;
  }

  PropMatrices::PrimaryMap 
  PropMatrices::LoadInterpMatrixMZ0(double M, double Z0) 
  {
    const double LgEmin = 12.;
    const double LgEmax = 22.;
    const double MaxDistance = 4.2122656250e+00 * utl::Gpc; 
    const int intM = 100*M;
    const int intZ0 = 100*Z0;
    std::string strM = std::to_string(abs(intM));
    std::string strZ0 = std::to_string(intZ0);

    TFile* fFile;
    std::string InterpFile = "./crp4/CRPropaG12_eEvoM";

    if( M < 0. ) 
      InterpFile += "m" + strM + "z" + strZ0 + "_nu.root";
    else 
      InterpFile += "p" + strM + "z" + strZ0 + "_nu.root";
  
    PropMatrices::PrimaryMap InterpM;
 
    fFile = TFile::Open(InterpFile.c_str());

    if (!fFile || fFile->IsZombie())
      throw runtime_error("cannot open file" + InterpFile);

    TIter nextkey = fFile->GetListOfKeys();
    TKey* key = 0;
    while ((key = dynamic_cast<TKey*>(nextkey()))) {
      TObject* const obj = key->ReadObj();
      if (obj->IsA()->InheritsFrom("TH2D")) {
        TH2D* h = static_cast<TH2D*>(obj);
        vector<string> splitname;
        stringstream  name(h->GetName());
        string line;
        while(getline(name, line, '_'))
	        splitname.push_back(line);
        if (splitname.size() != 3)
          throw runtime_error("cannot decode" + name.str());
        const unsigned int Aprim = stoi(splitname[1]);
        const unsigned int Asec = stoi(splitname[2]);
        TMatrixD& m = InterpM[Aprim][Asec];
        m.ResizeTo(h->GetNbinsY(), h->GetNbinsX());
        if (h->GetXaxis()->GetXmin() != LgEmin ) 
	        throw runtime_error("Different minimum energy, not supported");
        if (h->GetXaxis()->GetXmax() != LgEmax ) 
	        throw runtime_error("Different maximum energy, not supported");

        for (int iPrim = 0; iPrim < h->GetNbinsX(); ++iPrim) {
          const double dESource =
            pow(10, h->GetXaxis()->GetBinUpEdge(iPrim+1)) -
            pow(10, h->GetXaxis()->GetBinLowEdge(iPrim+1));
          for (int jSec = 0; jSec < h->GetNbinsY(); ++jSec) {
            const double dEEarth =
              pow(10, h->GetYaxis()->GetBinUpEdge(jSec+1)) -
              pow(10, h->GetYaxis()->GetBinLowEdge(jSec+1));
            m[jSec][iPrim] = h->GetBinContent(iPrim+1, jSec+1) * dESource / dEEarth;
          }
        }
      }
      else if (obj->IsA()->InheritsFrom("TH1D")) {
          TH1D* h = static_cast<TH1D*>(obj);
          if (string(h->GetName()) == string("hSimSettings")) {
            if( h->GetBinContent(1) != MaxDistance )
		          throw runtime_error("Different maximum distance, not supported");
          }
	    }
    }

    fFile->Close();
    fFile = nullptr;

    return InterpM;
  }

  void 
  PropMatrices::UpdateMZ0(double oldM, double newM, double oldZ0, double newZ0)
  {
    if( newM > Mmax || newM < Mmin )  
	    throw runtime_error("M drifted outside range: M = " + std::to_string(newM));
    if( newZ0 > Z0max || newZ0 < Z0min )
	    throw runtime_error("Z0 drifted outside range: Z0 = " + std::to_string(newZ0));

    // if nothing has changed don't update
    if(oldM == newM && oldZ0 == newZ0)
      return;

    const int DeltaM = floor(newM/dM)-floor(oldM/dM);
    const int DeltaZ0 = floor(newZ0/dZ0)-floor(oldZ0/dZ0); 
    const double xL = dM*floor(newM/dM);
    const double xR = xL + dM;
    const double yD = dZ0*floor(newZ0/dZ0);
    const double yU = yD + dZ0;

    // if in new region update all grid points
    if( abs(DeltaM) > 1 || abs(DeltaZ0) > 1 ) {
      InterpMatrixUL = PropMatrices::LoadInterpMatrixMZ0(xL, yU); 
      InterpMatrixUR = PropMatrices::LoadInterpMatrixMZ0(xR, yU); 
      InterpMatrixDL = PropMatrices::LoadInterpMatrixMZ0(xL, yD); 
      InterpMatrixDR = PropMatrices::LoadInterpMatrixMZ0(xR, yD); 
    }	
    // if new region shares corner with old region only bring in 3 new grid points
    else if( abs(DeltaM) == 1 && abs(DeltaZ0) == 1 ) {
      // ++
      if( DeltaM > 0 && DeltaZ0 > 0 ) {
        InterpMatrixDL = InterpMatrixUR;
        InterpMatrixUL = PropMatrices::LoadInterpMatrixMZ0(xL, yU); 
        InterpMatrixUR = PropMatrices::LoadInterpMatrixMZ0(xR, yU); 
        InterpMatrixDR = PropMatrices::LoadInterpMatrixMZ0(xR, yD); 
      }
      // +-
      else if( DeltaM > 0 && DeltaZ0 < 0 ) {
        InterpMatrixUL = InterpMatrixDR;
        InterpMatrixUR = PropMatrices::LoadInterpMatrixMZ0(xR, yU); 
        InterpMatrixDL = PropMatrices::LoadInterpMatrixMZ0(xL, yD); 
        InterpMatrixDR = PropMatrices::LoadInterpMatrixMZ0(xR, yD); 
      }
      // -+
      else if( DeltaM < 0 && DeltaZ0 > 0 ) {
        InterpMatrixDR = InterpMatrixUL;
        InterpMatrixUL = PropMatrices::LoadInterpMatrixMZ0(xL, yU); 
        InterpMatrixUR = PropMatrices::LoadInterpMatrixMZ0(xR, yU); 
        InterpMatrixDL = PropMatrices::LoadInterpMatrixMZ0(xL, yD); 
      }
      // --
      else {
        InterpMatrixUR = InterpMatrixDL;
        InterpMatrixUL = PropMatrices::LoadInterpMatrixMZ0(xL, yU); 
        InterpMatrixDL = PropMatrices::LoadInterpMatrixMZ0(xL, yD); 
        InterpMatrixDR = PropMatrices::LoadInterpMatrixMZ0(xR, yD); 
      }
    }
    // if new region is adjacent to old region only bring in 2 new grid points
    else if( (abs(DeltaM) == 1) != (abs(DeltaZ0) == 1) ) {
      // +M
      if( DeltaM == 1 ) {
        InterpMatrixUL = InterpMatrixUR;
        InterpMatrixDL = InterpMatrixDR;
        InterpMatrixUR = PropMatrices::LoadInterpMatrixMZ0(xR, yU); 
        InterpMatrixDR = PropMatrices::LoadInterpMatrixMZ0(xR, yD); 
      }
      // -M
      else if( DeltaM == -1 ) {
        InterpMatrixUR = InterpMatrixUL;
        InterpMatrixDR = InterpMatrixDL;
        InterpMatrixUL = PropMatrices::LoadInterpMatrixMZ0(xL, yU); 
        InterpMatrixDL = PropMatrices::LoadInterpMatrixMZ0(xL, yD); 
      }
      // +Z0
      else if( DeltaZ0 == 1 ) {
        InterpMatrixDL = InterpMatrixUL;
        InterpMatrixDR = InterpMatrixUR;
        InterpMatrixUL = PropMatrices::LoadInterpMatrixMZ0(xL, yU); 
        InterpMatrixUR = PropMatrices::LoadInterpMatrixMZ0(xR, yU); 
      }
      // -Z0
      else { 
        InterpMatrixUL = InterpMatrixDL;
        InterpMatrixUR = InterpMatrixDR;
        InterpMatrixDL = PropMatrices::LoadInterpMatrixMZ0(xL, yD); 
        InterpMatrixDR = PropMatrices::LoadInterpMatrixMZ0(xR, yD); 
      }
    }

    // if in same region then no grid points need updating

    // perform interpolation    
    const double x = newM;
    const double y = newZ0;

    InterpMZ0(x, xL, xR, y, yD, yU);

    return;
  }

  void 
  PropMatrices::InterpInitMZ0(double M, double Z0)
  {
    if( M > Mmax || M < Mmin )  
	    throw runtime_error("Initial M outside range: M = " + std::to_string(M));
    if( Z0 > Z0max || Z0 < Z0min )
	    throw runtime_error("Initial Z0 outside range: Z0 = " + std::to_string(Z0));

    const double x = M;
    const double y = Z0;
    const double xL = dM*floor(M/dM);
    const double xR = xL + dM;
    const double yD = dZ0*floor(Z0/dZ0);
    const double yU = yD + dZ0;

    InterpMatrixUL = PropMatrices::LoadInterpMatrixMZ0(xL, yU); 
    InterpMatrixUR = PropMatrices::LoadInterpMatrixMZ0(xR, yU); 
    InterpMatrixDL = PropMatrices::LoadInterpMatrixMZ0(xL, yD); 
    InterpMatrixDR = PropMatrices::LoadInterpMatrixMZ0(xR, yD); 

    SetEnergyRange(12., 22.);
    SetMaximumDistance(4.2122656250e+00 * utl::Gpc); 

    InterpMZ0(x, xL, xR, y, yD, yU);

    return;
  } 
  
  void
  PropMatrices::InterpMZ0(double x, double xL, double xR, double y, double yD, double yU)
  {
    ResetPrimaryMap();
    std::map<int, int> primkeyMap;
    std::map<int, int> seckeyMap;
    for ( auto const& iter1 : InterpMatrixUL) 
      primkeyMap[iter1.first] = 0;
    for ( auto const& iter1 : InterpMatrixUR) 
      primkeyMap[iter1.first] = 0; 
    for ( auto const& iter1 : InterpMatrixDL) 
      primkeyMap[iter1.first] = 0;
    for ( auto const& iter1 : InterpMatrixDR) 
      primkeyMap[iter1.first] = 0; 
    for ( auto const& iter1 : primkeyMap) {
      const unsigned int Aprim = iter1.first;
      seckeyMap.clear();
      for ( auto const& iter2 : InterpMatrixUL[Aprim]) 
        seckeyMap[iter2.first] = 0;
      for ( auto const& iter2 : InterpMatrixUR[Aprim]) 
        seckeyMap[iter2.first] = 0;
      for ( auto const& iter2 : InterpMatrixDL[Aprim]) 
        seckeyMap[iter2.first] = 0;
      for ( auto const& iter2 : InterpMatrixDR[Aprim]) 
        seckeyMap[iter2.first] = 0;
      for ( auto const& iter2 : seckeyMap) {
        const unsigned int Asec = iter2.first;
        TMatrixD& m = fMatrices[Aprim][Asec];
        const TMatrixD& mUL = InterpMatrixUL[Aprim][Asec];  
        const TMatrixD& mUR = InterpMatrixUR[Aprim][Asec];  
        const TMatrixD& mDL = InterpMatrixDL[Aprim][Asec];  
        const TMatrixD& mDR = InterpMatrixDR[Aprim][Asec];  
        const int maxNcols = max(mUL.GetNcols(), max(mUR.GetNcols(), max(mDL.GetNcols(), mDR.GetNcols())));
        
        if( (maxNcols != mUL.GetNcols() && mUL.GetNcols() > 0) || (maxNcols != mUR.GetNcols() && mUR.GetNcols() > 0) || 
            (maxNcols != mDL.GetNcols() && mDL.GetNcols() > 0) || (maxNcols != mDR.GetNcols()  && mDR.GetNcols() > 0) )
          throw runtime_error("PropMatrices column dimension mismatch!");
        
        const int maxNrows = max(mUL.GetNrows(), max(mUR.GetNrows(), max(mDL.GetNrows(), mDR.GetNrows())));
        
        if( (maxNrows != mUL.GetNrows() && mUL.GetNrows() > 0) || (maxNrows != mUR.GetNrows() && mUR.GetNrows() > 0) || 
            (maxNrows != mDL.GetNrows() && mDL.GetNrows() > 0) || (maxNrows != mDR.GetNrows()  && mDR.GetNrows() > 0) )
          throw runtime_error("PropMatrices row dimension mismatch!");
        
        m.ResizeTo(mUL.GetNrows(), mUL.GetNcols());
      //	fMatrices.SetEnergyRange(h->GetXaxis()->GetXmin(),
      //			       h->GetXaxis()->GetXmax());
        for (int iPrim = 0; iPrim < m.GetNcols(); ++iPrim) {
          for (int jSec = 0; jSec < m.GetNrows(); ++jSec) {
            const auto valUL = (mUL.GetNcols() == 0 || mUL.GetNrows() == 0)? 0. : mUL[jSec][iPrim];
            const auto valUR = (mUR.GetNcols() == 0 || mUR.GetNrows() == 0)? 0. : mUR[jSec][iPrim];
            const auto valDL = (mDL.GetNcols() == 0 || mDL.GetNrows() == 0)? 0. : mDL[jSec][iPrim];
            const auto valDR = (mDR.GetNcols() == 0 || mDR.GetNrows() == 0)? 0. : mDR[jSec][iPrim];
            // performs simple bicubic interpolation with zero derivative on edges of interpolation region (to ensure smooth derivatives at edges)
            const double X = (x-xL)/(xR-xL); 
            const double Y = (y-yD)/(yU-yD);
            const double a00 = valDL; 
            const double a22 = 9.*(valDL+valUR-valUL-valDR); 
            const double a33 = 4.*(valDL+valUR-valUL-valDR);
            const double a02 = 3.*(valUL-valDL); 
            const double a20 = 3.*(valDR-valDL);
            const double a03 = 2.*(valDL-valUL); 
            const double a30 = 2.*(valDL-valDR);
            const double a23 = 6.*(valDR+valUL-valDL-valUR); 
            const double a32 = a23;
            m[jSec][iPrim] = a00 + a02*pow(Y, 2) + a20*pow(X, 2) + a03*pow(Y, 3) + a30*pow(X, 3) +
                              a22*pow(X*Y, 2) + a33*pow(X*Y, 3) + a23*pow(X, 2)*pow(Y, 3) + a32*pow(X, 3)*pow(Y, 2);
           }
        }
       }
     }

    return;
  }

  PropMatrices::PrimaryMap 
  PropMatrices::LoadInterpMatrixDmin(double Dmin) 
  {
    const double LgEmin = 12.;
    const double LgEmax = 22.;
    const double MaxDistance = 4.2122656250e+00 * utl::Gpc; 
    std::string strDmin = std::to_string(int(Dmin));

    TFile* fFile;
    std::string InterpFile = "./Data/crp5/CRPropaG12_SFR2_";
    InterpFile += strDmin + "_nu.root";
  
    PropMatrices::PrimaryMap InterpM;
 
    fFile = TFile::Open(InterpFile.c_str());

    if (!fFile || fFile->IsZombie())
      throw runtime_error("cannot open file" + InterpFile);

    TIter nextkey = fFile->GetListOfKeys();
    TKey* key = 0;
    while ((key = dynamic_cast<TKey*>(nextkey()))) {
      TObject* const obj = key->ReadObj();
      if (obj->IsA()->InheritsFrom("TH2D")) {
        TH2D* h = static_cast<TH2D*>(obj);
        vector<string> splitname;
        stringstream  name(h->GetName());
        string line;
        while(getline(name, line, '_'))
	        splitname.push_back(line);
        if (splitname.size() != 3)
          throw runtime_error("cannot decode" + name.str());
        const unsigned int Aprim = stoi(splitname[1]);
        const unsigned int Asec = stoi(splitname[2]);
        TMatrixD& m = InterpM[Aprim][Asec];
        m.ResizeTo(h->GetNbinsY(), h->GetNbinsX());
        if (h->GetXaxis()->GetXmin() != LgEmin ) 
	        throw runtime_error("Different minimum energy, not supported");
        if (h->GetXaxis()->GetXmax() != LgEmax ) 
	        throw runtime_error("Different maximum energy, not supported");

        for (int iPrim = 0; iPrim < h->GetNbinsX(); ++iPrim) {
          const double dESource =
            pow(10, h->GetXaxis()->GetBinUpEdge(iPrim+1)) -
            pow(10, h->GetXaxis()->GetBinLowEdge(iPrim+1));
          for (int jSec = 0; jSec < h->GetNbinsY(); ++jSec) {
            const double dEEarth =
              pow(10, h->GetYaxis()->GetBinUpEdge(jSec+1)) -
              pow(10, h->GetYaxis()->GetBinLowEdge(jSec+1));
            m[jSec][iPrim] = h->GetBinContent(iPrim+1, jSec+1) * dESource / dEEarth;
          }
        }
      }
      else if (obj->IsA()->InheritsFrom("TH1D")) {
          TH1D* h = static_cast<TH1D*>(obj);
          if (string(h->GetName()) == string("hSimSettings")) {
            if( h->GetBinContent(1) != MaxDistance )
		          throw runtime_error("Different maximum distance, not supported");
          }
	    }
    }

    fFile->Close();
    fFile = nullptr;

    return InterpM;
  }

  void 
  PropMatrices::UpdateDmin(double oldDmin, double newDmin)
  {
    if( newDmin > DminMax || newDmin < DminMin )  
	    throw runtime_error("Dmin drifted outside range: Dmin = " + std::to_string(newDmin) + " Mpc");

    // if nothing has changed don't update
    if(oldDmin == newDmin) 
      return;

    const int newposR = (upper_bound(DminGrid.begin(), DminGrid.end(), newDmin) - DminGrid.begin());
    const int deltapos = newposR - posR;
    const double xL = DminGrid[newposR-1];
    const double xR = DminGrid[newposR];

    // if in new region update all grid points
    if( abs(deltapos) > 1 ) {
      InterpMatrixL = PropMatrices::LoadInterpMatrixDmin(xL); 
      InterpMatrixR = PropMatrices::LoadInterpMatrixDmin(xR); 
    }	
    // if in adjacent region update one grid point 
    else if( abs(deltapos) == 1 ) {
      // +
      if( deltapos > 0 ) {
        InterpMatrixL = InterpMatrixR;
        InterpMatrixR = PropMatrices::LoadInterpMatrixDmin(xR); 
      }
      // -
      else {
	      InterpMatrixR = InterpMatrixL;
	      InterpMatrixL = PropMatrices::LoadInterpMatrixDmin(xL); 
      }
    }
    // if in same region just update interpolation 

    // perform interpolation    

    const double x = newDmin;
    InterpDmin(x, xL, xR);
    posR = newposR;

    return;
  }

  void 
  PropMatrices::InterpInitDmin(double Dmin) 
  {
    if( Dmin > DminMax || Dmin < DminMin )  
	    throw runtime_error("Initial Dmin outside range: Dmin = " + std::to_string(Dmin) + " Mpc");

    posR = (upper_bound(DminGrid.begin(), DminGrid.end(), Dmin) - DminGrid.begin());
    const double x = Dmin;
    const double xL = DminGrid[posR-1];
    const double xR = DminGrid[posR];

    InterpMatrixL = PropMatrices::LoadInterpMatrixDmin(xL); 
    InterpMatrixR = PropMatrices::LoadInterpMatrixDmin(xR); 

    SetEnergyRange(12., 22.);
    SetMaximumDistance(4.2122656250e+00 * utl::Gpc); 
    
    InterpDmin(x, xL, xR);

    return;
  } 

  void 
  PropMatrices::InterpDmin(double x, double xL, double xR)
  {
    ResetPrimaryMap();
    std::map<int, int> primkeyMap;
    std::map<int, int> seckeyMap;
    for ( auto const& iter1 : InterpMatrixL) 
      primkeyMap[iter1.first] = 0;
    for ( auto const& iter1 : InterpMatrixR) 
      primkeyMap[iter1.first] = 0; 
    for ( auto const& iter1 : primkeyMap) {
      const unsigned int Aprim = iter1.first;
      seckeyMap.clear();
      for ( auto const& iter2 : InterpMatrixL[Aprim]) 
        seckeyMap[iter2.first] = 0;
      for ( auto const& iter2 : InterpMatrixR[Aprim]) 
        seckeyMap[iter2.first] = 0;
      for ( auto const& iter2 : seckeyMap) {
        const unsigned int Asec = iter2.first;
        TMatrixD& m = GetMatrix(Aprim, Asec);
        const TMatrixD& mL = InterpMatrixL[Aprim][Asec];  
        const TMatrixD& mR = InterpMatrixR[Aprim][Asec];
        if( mL.GetNcols() != mR.GetNcols() && mL.GetNcols() > 0 && mR.GetNcols() > 0 )
          throw runtime_error("PropMatrices column dimension mismatch!");
        if( mL.GetNrows() != mR.GetNrows() && mL.GetNrows() > 0 && mR.GetNrows() > 0 )
          throw runtime_error("PropMatrices row dimension mismatch!");
        m.ResizeTo(mL.GetNrows(), mL.GetNcols());
        for (int iPrim = 0; iPrim < m.GetNcols(); ++iPrim) {
          for (int jSec = 0; jSec < m.GetNrows(); ++jSec) {
            const auto valL = (mL.GetNcols() == 0 || mL.GetNrows() == 0)? 0. : mL[jSec][iPrim];
            const auto valR = (mR.GetNcols() == 0 || mR.GetNrows() == 0)? 0. : mR[jSec][iPrim];
            // performs simple cubic interpolation with zero derivative at endpoints (to ensure smooth derivatives at gridpoints)
            const double X = (x-xL)/(xR-xL);
            const double C1 = -2.*(valR-valL);
            const double C2 = 3.*(valR-valL);
            const double C3 = valL;
            m[jSec][iPrim] = C1*pow(X, 3) + C2*pow(X, 2) + C3;
          }
        }
      }
    }

    return;
  }
}
