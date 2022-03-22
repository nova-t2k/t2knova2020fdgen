#include "ChannelHistCollections.h"

#include "T2KNOvA/TrueSelectionHelper.hxx"

#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TTreeReaderValue.h"

#include "TLorentzVector.h"

namespace t2knova {

struct vect4 {
  vect4(std::array<float, 4> const &vec) {
    TLorentzVector v(vec[0], vec[1], vec[2], vec[3]);
    E = v.E();
    _gamma = v.Gamma();
  }
  float E;
  float _gamma;
  float Gamma() { return _gamma; }
};

struct nova_part {
  nova_part(int _pdg, std::array<float, 4> const &vec) : pdg(_pdg), p(vec) {}
  int pdg;
  vect4 p;
};

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

  float Eav_NOvA() {

    float eAvail = 0;

    for (int i = 0; i < *_nfsp; ++i) {
      nova_part particle(_pdg[i],
                         std::array<float, 4>{_px[i], _py[i], _pz[i], _E[i]});

      // Code below here is from Zoya (unaltered from NOvA/CAFAna
      // implementation)

      double particle_Eavail = 0;
      // protons, pions: KE only
      if (particle.pdg == 2212 || abs(particle.pdg) == 211) {
        double gamma = particle.p.Gamma();
        particle_Eavail = (gamma - 1) / gamma * particle.p.E;
      }
      // pi0s, electrons, photons: total energy
      else if (particle.pdg == 111 || abs(particle.pdg) == 11 ||
               particle.pdg == 22)
        particle_Eavail = particle.p.E;

      else if (particle.pdg >= 2000000000) {
        // skip the bindinos
      }

      else if (particle.pdg >= 1000000000) {
        // do nothing for nucleons
      } else if (particle.pdg >= 2000 && particle.pdg != 2212 &&
                 particle.pdg != 2112) {
        particle_Eavail = particle.p.E - 0.9382;
        // Primarily strange baryons add total energy minus proton mass since
        // decays will mostly contain protons
      } else if (particle.pdg <= -2000) {
        particle_Eavail = particle.p.E + 0.9382;
        // Primarily anti-protons add total energy plus proton mass since
        // anhillation is mostly the interaction mode
      } else if (particle.pdg != 2112 &&
                 (abs(particle.pdg) < 11 ||
                  abs(particle.pdg) > 16)) { // no neutrons or leptons
        particle_Eavail = particle.p.E;
        // mostly kaons add all the energy
      }

      eAvail += particle_Eavail;
    } // end loop over primaries

    return eAvail;
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

  std::vector<int> GetSelections(int Mode = 0) {
    return t2knova::GetSelections(
        T2KNOvAFlatTreeToFSParticleSummary(*_nfsp, (int *)_pdg.GetAddress(),
                                           (float *)_E.GetAddress()),
        Mode);
  }

  int GetPrimarySelection(int Mode = 0) {
    return t2knova::GetPrimarySelection(
        T2KNOvAFlatTreeToFSParticleSummary(*_nfsp, (int *)_pdg.GetAddress(),
                                           (float *)_E.GetAddress()),
        Mode);
  }

private:
  // I know, we all hate it, but because of design of TTreeReader I really think
  // having getters is easier
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