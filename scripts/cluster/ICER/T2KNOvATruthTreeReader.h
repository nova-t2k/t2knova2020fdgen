#include "ChannelHistCollections.h"
#include "T2KNOvATrueSelectionHelper.hxx"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TTreeReaderValue.h"

namespace t2knova {

class T2KNOvATruthTreeReader {
 public:
  int PDGLep() { return *_PDGLep; }
  int PDGNu() { return *_PDGNu; }
  int Mode() { return fModeInfo ? *_Mode : 0; }
  float CosLep() { return *_CosLep; }
  float EavAlt() { return *_EavAlt; }
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

  void SetNoModeInfo() { fModeInfo = false; }

  T2KNOvATruthTreeReader(TTreeReader& rdr)
      : _PDGLep(rdr, "PDGLep"),
        _PDGNu(rdr, "PDGnu"),
        _Mode(rdr, "Mode"),
        _CosLep(rdr, "CosLep"),
        _EavAlt(rdr, "EavAlt"),
        _Enu_true(rdr, "Enu_true"),
        _PLep(rdr, "PLep"),
        _Q2(rdr, "Q2"),
        _q0(rdr, "q0"),
        _q3(rdr, "q3"),
        _tgta(rdr, "tgta"),
        _nfsp(rdr, "nfsp"),
        _pdg(rdr, "pdg"),
        _E(rdr, "E"),
        _fScaleFactor(rdr, "fScaleFactor"),
        _RWWeight(rdr, "RWWeight"),
        fModeInfo(true) {
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
    CheckValue(_E);
    CheckValue(_fScaleFactor);
    CheckValue(_RWWeight);
    rdr.Restart();
  }

  std::vector<int> GetSelections() {
    return t2knova::GetSelections(
        T2KNOvAFlatTreeToFSParticleSummary(*_nfsp, (int*)_pdg.GetAddress()));
  }

  int GetPrimarySelection() {
    return t2knova::GetPrimarySelection(
        T2KNOvAFlatTreeToFSParticleSummary(*_nfsp, (int*)_pdg.GetAddress()));
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

  TTreeReaderValue<double> _fScaleFactor;
  TTreeReaderValue<double> _RWWeight;

  bool fModeInfo;

  bool CheckValue(ROOT::Internal::TTreeReaderValueBase& value) {
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

}  // namespace t2knova