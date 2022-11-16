#include "ChannelHistCollections.h"
#include "NOvAFuncs.h"

#include "T2KNOvA/TrueSelectionHelper.hxx"

#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TTreeReaderValue.h"

#include "TLorentzVector.h"

namespace t2knova {

class T2KNOvATruthTreeReader {
public:
  int PDGLep() { return *_PDGLep; }
  int PDGNu() { return *_PDGNu; }
  int Mode() { return fModeInfo ? *_Mode : 0; }
  float CosLep() { return *_CosLep; }
  float AngLep_deg() { return acos(*_CosLep) * (180.0 / M_PI); }
  float EavAlt() { return *_EavAlt; }

  TLorentzVector FSLepP4() {
    for (int i = 0; i < *_nfsp; ++i) {
      if (_pdg[i] == PDGLep()) {
        return TLorentzVector(_px[i], _py[i], _pz[i], _E[i]);
      }
    }
    return TLorentzVector();
  }

  float hmfscpip() {
    double maxp = 0;
    for (int i = 0; i < *_nfsp; ++i) {
      if (std::abs(_pdg[i]) == 211) {
        TVector3 v(_px[i], _py[i], _pz[i]);
        if (v.Mag() > maxp) {
          maxp = v.Mag();
        }
      }
    }
    return maxp;
  }
  float hmfspi0p() {
    double maxp = 0;
    for (int i = 0; i < *_nfsp; ++i) {
      if (std::abs(_pdg[i]) == 111) {
        TVector3 v(_px[i], _py[i], _pz[i]);
        if (v.Mag() > maxp) {
          maxp = v.Mag();
        }
      }
    }
    return maxp;
  }
  int ncpi() {
    int ncpi = 0;
    for (int i = 0; i < *_nfsp; ++i) {
      if (std::abs(_pdg[i]) == 211) {
        ncpi++;
      }
    }
    return ncpi;
  }
  int npi0() {
    int npi0 = 0;
    for (int i = 0; i < *_nfsp; ++i) {
      if (std::abs(_pdg[i]) == 111) {
        npi0++;
      }
    }
    return npi0;
  }

  float hmfsprotonp() {
    double maxp = 0;
    for (int i = 0; i < *_nfsp; ++i) {
      if (_pdg[i] == 2212) {
        TVector3 v(_px[i], _py[i], _pz[i]);
        if (v.Mag() > maxp) {
          maxp = v.Mag();
        }
      }
    }
    return maxp;
  }
  float hmfsneutronp() {
    double maxp = 0;
    for (int i = 0; i < *_nfsp; ++i) {
      if (_pdg[i] == 2112) {
        TVector3 v(_px[i], _py[i], _pz[i]);
        if (v.Mag() > maxp) {
          maxp = v.Mag();
        }
      }
    }
    return maxp;
  }
  int nproton() {
    int nproton = 0;
    for (int i = 0; i < *_nfsp; ++i) {
      if (_pdg[i] == 2212) {
        nproton++;
      }
    }
    return nproton;
  }
  int nneutron() {
    int nneutron = 0;
    for (int i = 0; i < *_nfsp; ++i) {
      if (_pdg[i] == 2112) {
        nneutron++;
      }
    }
    return nneutron;
  }

  float Eav_NOvA() {
    return t2knova::Eav_NOvA(*_nfsp, &_pdg[0], &_px[0], &_py[0], &_pz[0],
                             &_E[0]);
  }

  selection NOvAFSIMode() {
    return t2knova::NOvAFSIMode(*_nfsp, &_pdg[0], &_px[0], &_py[0], &_pz[0],
                                &_E[0]);
  }

  float Enu_true() { return *_Enu_true; }
  float PLep() { return *_PLep; }
  float Q2() { return *_Q2; }
  float q0() { return *_q0; }
  float q3() { return *_q3; }
  int tgta() { return *_tgta; }
  float fScaleFactor() { return *_fScaleFactor; }
  float RWWeight() { return *_RWWeight; }

  float EGamma() {
    float Esum = 0;
    for (int i = 0; i < *_nfsp; ++i) {
      if (_pdg[i] == 22) {
        Esum += _E[i];
      }
    }
    return Esum;
  }

  // BANFFEventBase::
  float GetTrueEnuQE() {
    static double const proton_mass = 0.938272;  // GeV
    static double const neutron_mass = 0.939566; // GeV
    bool nu = (PDGNu() > 0);
    double Eb = 27. * 0.001; // GeV
    double m_in = neutron_mass - Eb * 1.e-3;
    double m_out = proton_mass;
    if (!nu) {
      m_in = proton_mass - Eb * 1.e-3;
      m_out = neutron_mass;
    }
    auto pfslep = FSLepP4();
    double muon_energy = pfslep.E();
    return (2. * m_in * muon_energy - pfslep.Mag2() - m_in * m_in +
            m_out * m_out) /
           (2. * (m_in - muon_energy +
                  pfslep.Vect().Mag() * pfslep.Vect().CosTheta()));
  }

  // BANFFEventBase::
  float GetQ2QE() {
    auto pfslep = FSLepP4();
    double muon_energy = pfslep.E();
    return -pfslep.Mag2() +
           2. * GetTrueEnuQE() *
               (muon_energy - pfslep.Vect().Mag() * pfslep.Vect().CosTheta());
  }

  std::string PrintStack() {
    std::stringstream ss;
    ss << "NFSP: " << (*_nfsp) << "\n";

    for (int i = 0; i < (*_nfsp); ++i) {
      ss << "\t" << i << ", PDG: " << _pdg[i] << ", E: " << _E[i] << "\n";
    }
    return ss.str();
  }

  void SetNoModeInfo() { fModeInfo = false; }

  T2KNOvATruthTreeReader(TTreeReader &rdr)
      : _PDGLep(rdr, "PDGLep"), _PDGNu(rdr, "PDGnu"), _Mode(rdr, "Mode"),
        _CosLep(rdr, "CosLep"), _EavAlt(rdr, "EavAlt"),
        _Enu_true(rdr, "Enu_true"), _PLep(rdr, "PLep"), _Q2(rdr, "Q2"),
        _q0(rdr, "q0"), _q3(rdr, "q3"), _tgta(rdr, "tgta"), _nfsp(rdr, "nfsp"),
        _pdg(rdr, "pdg"), _px(rdr, "px"), _py(rdr, "py"), _pz(rdr, "pz"),
        _E(rdr, "E"), _fScaleFactor(rdr, "fScaleFactor"),
        _RWWeight(rdr, "RWWeight"), fModeInfo(true) {
    rdr.Restart();
    (void)rdr.Next();
    CheckValue(_PDGLep);
    CheckValue(_PDGNu);
    CheckValue(_Mode);
    CheckValue(_CosLep);
    CheckValue(_EavAlt);
    CheckValue(_Enu_true);
    CheckValue(_PLep);
    CheckValue(_Q2);
    CheckValue(_q0);
    CheckValue(_q3);
    CheckValue(_tgta);
    CheckValue(_nfsp);
    CheckValue(_pdg);
    CheckValue(_px);
    CheckValue(_py);
    CheckValue(_pz);
    CheckValue(_E);
    CheckValue(_fScaleFactor);
    CheckValue(_RWWeight);
    rdr.Restart();
  }

  std::vector<int> GetSelections() {
    return t2knova::GetSelections(T2KNOvAFlatTreeToFSParticleSummary(
                                      *_nfsp, (int *)_pdg.GetAddress(),
                                      (float *)_E.GetAddress(), NOvAFSIMode()),
                                  *_Mode);
  }

  int GetPrimarySelection() {
    return t2knova::GetPrimarySelection(
        T2KNOvAFlatTreeToFSParticleSummary(*_nfsp, (int *)_pdg.GetAddress(),
                                           (float *)_E.GetAddress(),
                                           NOvAFSIMode()),
        *_Mode);
  }

private:
  // I know, we all hate it, but because of design of TTreeReader I really
  // think having getters is easier
  TTreeReaderValue<int> _PDGLep;
  TTreeReaderValue<int> _PDGNu;
  TTreeReaderValue<int> _Mode;

  TTreeReaderValue<float> _CosLep;
  TTreeReaderValue<float> _EavAlt;
  TTreeReaderValue<float> _Enu_true;
  TTreeReaderValue<float> _PLep;
  TTreeReaderValue<float> _Q2;
  TTreeReaderValue<float> _q0;
  TTreeReaderValue<float> _q3;
  TTreeReaderValue<int> _tgta;

  TTreeReaderValue<int> _nfsp;
  TTreeReaderArray<int> _pdg;
  TTreeReaderArray<float> _E;
  TTreeReaderArray<float> _px;
  TTreeReaderArray<float> _py;
  TTreeReaderArray<float> _pz;

  TTreeReaderValue<double> _fScaleFactor;
  TTreeReaderValue<double> _RWWeight;

  bool fModeInfo;

  bool CheckValue(ROOT::Internal::TTreeReaderValueBase &value) {
    if (value.GetSetupStatus() < 0) {
      std::cerr << "Error " << value.GetSetupStatus()
                << "setting up reader for " << value.GetBranchName() << '\n';
      return false;
    }
    // else {
    //   std::cout << "Successfully read branch: " << value.GetBranchName()
    //             << std::endl;
    // }
    return true;
  }
};

} // namespace t2knova
