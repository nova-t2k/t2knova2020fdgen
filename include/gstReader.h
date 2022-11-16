#include "ChannelHistCollections.h"
#include "NOvAFuncs.h"

#include "T2KNOvA/TrueSelectionHelper.hxx"

#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TTreeReaderValue.h"

#include "TLorentzVector.h"

namespace t2knova {

class gstReader {
public:
  int PDGLep() { return *_PDGLep; }
  int PDGNu() { return *_PDGNu; }
  int Mode() { return fModeInfo ? *_Mode : 0; }

  TLorentzVector FSLepP4() {
    FillCompatInfo();

    for (int i = 0; i < _compat_nfsp; ++i) {
      if (_compat_pdg[i] == PDGLep()) {
        return TLorentzVector(_compat_px[i], _compat_py[i], _compat_pz[i],
                              _compat_E[i]);
      }
    }
    return TLorentzVector();
  }

  float hmfscpip() {
    FillCompatInfo();

    double maxp = 0;
    for (int i = 0; i < _compat_nfsp; ++i) {
      if (std::abs(_compat_pdg[i]) == 211) {
        TVector3 v(_compat_px[i], _compat_py[i], _compat_pz[i]);
        if (v.Mag() > maxp) {
          maxp = v.Mag();
        }
      }
    }
    return maxp;
  }
  float hmfspi0p() {
    FillCompatInfo();

    double maxp = 0;
    for (int i = 0; i < _compat_nfsp; ++i) {
      if (std::abs(_compat_pdg[i]) == 111) {
        TVector3 v(_compat_px[i], _compat_py[i], _compat_pz[i]);
        if (v.Mag() > maxp) {
          maxp = v.Mag();
        }
      }
    }
    return maxp;
  }
  int ncpi() {
    FillCompatInfo();

    int ncpi = 0;
    for (int i = 0; i < _compat_nfsp; ++i) {
      if (std::abs(_compat_pdg[i]) == 211) {
        ncpi++;
      }
    }
    return ncpi;
  }
  int npi0() {
    FillCompatInfo();

    int npi0 = 0;
    for (int i = 0; i < _compat_nfsp; ++i) {
      if (std::abs(_compat_pdg[i]) == 111) {
        npi0++;
      }
    }
    return npi0;
  }

  float hmfsprotonp() {
    FillCompatInfo();

    double maxp = 0;
    for (int i = 0; i < _compat_nfsp; ++i) {
      if (_compat_pdg[i] == 2212) {
        TVector3 v(_compat_px[i], _compat_py[i], _compat_pz[i]);
        if (v.Mag() > maxp) {
          maxp = v.Mag();
        }
      }
    }
    return maxp;
  }
  float hmfsneutronp() {
    FillCompatInfo();

    double maxp = 0;
    for (int i = 0; i < _compat_nfsp; ++i) {
      if (_compat_pdg[i] == 2112) {
        TVector3 v(_compat_px[i], _compat_py[i], _compat_pz[i]);
        if (v.Mag() > maxp) {
          maxp = v.Mag();
        }
      }
    }
    return maxp;
  }
  int nproton() {
    FillCompatInfo();

    int nproton = 0;
    for (int i = 0; i < _compat_nfsp; ++i) {
      if (_compat_pdg[i] == 2212) {
        nproton++;
      }
    }
    return nproton;
  }
  int nneutron() {
    FillCompatInfo();

    int nneutron = 0;
    for (int i = 0; i < _compat_nfsp; ++i) {
      if (_compat_pdg[i] == 2112) {
        nneutron++;
      }
    }
    return nneutron;
  }

  float Eav_NOvA() {
    FillCompatInfo();

    return t2knova::Eav_NOvA(_compat_nfsp, _compat_pdg.data(),
                             _compat_px.data(), _compat_py.data(),
                             _compat_pz.data(), _compat_E.data());
  }

  selection NOvAFSIMode() {
    FillCompatInfo();

    return t2knova::NOvAFSIMode(_compat_nfsp, _compat_pdg.data(),
                                _compat_px.data(), _compat_py.data(),
                                _compat_pz.data(), _compat_E.data());
  }

  float Enu_true() { return *_Enu_true; }
  float PLep() { return *_PLep; }
  int tgta() { return *_tgta; }

  float EGamma() {
    FillCompatInfo();

    float Esum = 0;
    for (int i = 0; i < _compat_nfsp; ++i) {
      if (_compat_pdg[i] == 22) {
        Esum += _compat_E[i];
      }
    }
    return Esum;
  }

  // BANFFEventBase::
  float GetTrueEnuQE() {
    FillCompatInfo();

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
    FillCompatInfo();

    auto pfslep = FSLepP4();
    double muon_energy = pfslep.E();
    return -pfslep.Mag2() +
           2. * GetTrueEnuQE() *
               (muon_energy - pfslep.Vect().Mag() * pfslep.Vect().CosTheta());
  }

  std::string PrintStack() {
    FillCompatInfo();

    std::stringstream ss;
    ss << "NFSP: " << (_compat_nfsp) << "\n";

    for (int i = 0; i < (_compat_nfsp); ++i) {
      ss << "\t" << i << ", PDG: " << _compat_pdg[i] << ", E: " << _compat_E[i]
         << "\n";
    }
    return ss.str();
  }

  void SetNoModeInfo() { fModeInfo = false; }

  gstReader(TTreeReader &rdr)
      : _iev(rdr, "iev"), compat_i(-1), _PDGLep(rdr, "fspl"), _ELep(rdr, "El"),
        _PXLep(rdr, "pxl"), _PYLep(rdr, "pyl"), _PZLep(rdr, "pzl"),
        _PDGNu(rdr, "neu"), _Mode(rdr, "neut_code"), _Enu_true(rdr, "Ev"),
        _PLep(rdr, "pl"), _tgta(rdr, "A"), _nfsp(rdr, "nf"), _pdg(rdr, "pdgf"),
        _px(rdr, "pxf"), _py(rdr, "pyf"), _pz(rdr, "pzf"), _E(rdr, "Ef"),
        fModeInfo(true) {
    rdr.Restart();
    (void)rdr.Next();
    CheckValue(_PDGLep);
    CheckValue(_PDGNu);
    CheckValue(_Mode);
    CheckValue(_Enu_true);
    CheckValue(_PLep);
    CheckValue(_tgta);
    CheckValue(_nfsp);
    CheckValue(_pdg);
    CheckValue(_px);
    CheckValue(_py);
    CheckValue(_pz);
    CheckValue(_E);
    rdr.Restart();
  }

  std::vector<int> GetSelections() {
    FillCompatInfo();

    return t2knova::GetSelections(
        T2KNOvAFlatTreeToFSParticleSummary(_compat_nfsp, _compat_pdg.data(),
                                           _compat_E.data(), NOvAFSIMode()),
        *_Mode);
  }

  int GetPrimarySelection() {
    FillCompatInfo();

    return t2knova::GetPrimarySelection(
        T2KNOvAFlatTreeToFSParticleSummary(_compat_nfsp, _compat_pdg.data(),
                                           _compat_E.data(), NOvAFSIMode()),
        *_Mode);
  }

private:
  // I know, we all hate it, but because of design of TTreeReader I really
  // think having getters is easier
  TTreeReaderValue<int> _iev;
  int compat_i;

  int _compat_nfsp;

  TTreeReaderValue<int> _PDGLep;
  TTreeReaderValue<double> _ELep;
  TTreeReaderValue<double> _PXLep;
  TTreeReaderValue<double> _PYLep;
  TTreeReaderValue<double> _PZLep;
  TTreeReaderValue<int> _PDGNu;
  TTreeReaderValue<int> _Mode;

  TTreeReaderValue<double> _Enu_true;
  TTreeReaderValue<double> _PLep;
  TTreeReaderValue<int> _tgta;

  TTreeReaderValue<int> _nfsp;
  TTreeReaderArray<int> _pdg;
  TTreeReaderArray<double> _E;
  TTreeReaderArray<double> _px;
  TTreeReaderArray<double> _py;
  TTreeReaderArray<double> _pz;

  std::vector<int> _compat_pdg;
  std::vector<double> _compat_E;
  std::vector<double> _compat_px;
  std::vector<double> _compat_py;
  std::vector<double> _compat_pz;

  void FillCompatInfo() {
    if (*_iev == compat_i) {
      return;
    }

    _compat_nfsp = *_nfsp + 1;
    _compat_pdg.clear();
    _compat_E.clear();
    _compat_px.clear();
    _compat_py.clear();
    _compat_pz.clear();

    _compat_pdg.push_back(*_PDGLep);
    _compat_E.push_back(*_ELep);
    _compat_px.push_back(*_PXLep);
    _compat_py.push_back(*_PYLep);
    _compat_pz.push_back(*_PZLep);

    for (int i = 0; i < *_nfsp; ++i) {
      _compat_pdg.push_back(_pdg[i]);
      _compat_E.push_back(_E[i]);
      _compat_px.push_back(_px[i]);
      _compat_py.push_back(_py[i]);
      _compat_pz.push_back(_pz[i]);
    }

    compat_i = *_iev;
  }

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
