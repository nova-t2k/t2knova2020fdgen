#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "THn.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TProfile.h"
#include "TStyle.h"
#include "TTreeReader.h"

struct FlagBlob {
  bool flagCCINC;
  bool flagCC0pi;
  bool flagCC1cpi;
  bool flagCC1pi0;
  bool flagNCINC;
  bool flagNC0pi;
  bool flagNC1cpi;
  bool flagNC1pi0;
};

template <typename TH> struct THTraits {};

template <> struct THTraits<TH1D> { using Base = TH1; };
template <> struct THTraits<TH2D> { using Base = TH1; };
template <> struct THTraits<TH3D> { using Base = TH1; };

template <> struct THTraits<THnD> { using Base = THnBase; };

template <typename TH> struct hblob {
  TH CCInc;
  TH CC0Pi;
  TH CC1CPi;
  TH CC1Pi0;
  TH CCOther;
  TH NCInc;
  TH NC0Pi;
  TH NC1CPi;
  TH NC1Pi0;
  TH NCOther;

  using THBase = typename THTraits<TH>::Base;

  hblob() {}

  template <typename... Args>
  hblob(std::string const &name, std::string const &title, Args... binning)
      : CCInc((name + "_CCInc").c_str(), title.c_str(), binning...),
        CC0Pi((name + "_CC0Pi").c_str(), title.c_str(), binning...),
        CC1CPi((name + "_CC1CPi").c_str(), title.c_str(), binning...),
        CC1Pi0((name + "_CC1Pi0").c_str(), title.c_str(), binning...),
        CCOther((name + "_CCOther").c_str(), title.c_str(), binning...),
        NCInc((name + "_NCInc").c_str(), title.c_str(), binning...),
        NC0Pi((name + "_NC0Pi").c_str(), title.c_str(), binning...),
        NC1CPi((name + "_NC1CPi").c_str(), title.c_str(), binning...),
        NC1Pi0((name + "_NC1Pi0").c_str(), title.c_str(), binning...),
        NCOther((name + "_NCOther").c_str(), title.c_str(), binning...) {

    Apply([=](TH &h) { h.SetDirectory(nullptr); });
    Apply([=](TH &h) { h.Sumw2(true); });
  }

  void SetName(std::string const &name) {
    CCInc.SetName((name + "_CCInc").c_str());
    CC0Pi.SetName((name + "_CC0Pi").c_str());
    CC1CPi.SetName((name + "_CC1CPi").c_str());
    CC1Pi0.SetName((name + "_CC1Pi0").c_str());
    CCOther.SetName((name + "_CCOther").c_str());
    NCInc.SetName((name + "_NCInc").c_str());
    NC0Pi.SetName((name + "_NC0Pi").c_str());
    NC1CPi.SetName((name + "_NC1CPi").c_str());
    NC1Pi0.SetName((name + "_NC1Pi0").c_str());
    NCOther.SetName((name + "_NCOther").c_str());
  }

  void SetTitle(std::string const &title) {
    Apply([=](TH &h) { h.SetTitle(title.c_str()); });
  }

  void Write(TDirectory *f, bool scale = false) {
    Apply([=](TH &h) {
      if (scale) {
        THBase *hc = (THBase *)h.Clone();
        hc->SetDirectory(nullptr);
        hc->Scale(1, "width");
        f->WriteTObject(hc, h.GetName());
        delete hc;
      } else {
        f->WriteTObject(&h, h.GetName());
      }
    });
  }

  template <typename... XY> void Fill(double w, FlagBlob const &blb, XY... xy) {
    if (blb.flagCCINC) {

      CCInc.Fill(xy..., w);

      if (blb.flagCC0pi) {
        CC0Pi.Fill(xy..., w);
      } else if (blb.flagCC1cpi) {
        CC1CPi.Fill(xy..., w);
      } else if (blb.flagCC1pi0) {
        CC1Pi0.Fill(xy..., w);
      } else {
        CCOther.Fill(xy..., w);
      }
    } else if (blb.flagNCINC) {

      NCInc.Fill(xy..., w);

      if (blb.flagNC0pi) {
        NC0Pi.Fill(xy..., w);
      } else if (blb.flagNC1cpi) {
        NC1CPi.Fill(xy..., w);
      } else if (blb.flagNC1pi0) {
        NC1Pi0.Fill(xy..., w);
      } else {
        NCOther.Fill(xy..., w);
      }
    }
  }

  template <typename THOut>
  hblob<THOut> Transform(std::function<THOut(TH const &)> f) {
    hblob<THOut> out;
    out.CCInc = THOut(f(CCInc));
    out.CC0Pi = THOut(f(CC0Pi));
    out.CC1CPi = THOut(f(CC1CPi));
    out.CC1Pi0 = THOut(f(CC1Pi0));
    out.CCOther = THOut(f(CCOther));
    out.NCInc = THOut(f(NCInc));
    out.NC0Pi = THOut(f(NC0Pi));
    out.NC1CPi = THOut(f(NC1CPi));
    out.NC1Pi0 = THOut(f(NC1Pi0));
    out.NCOther = THOut(f(NCOther));
    return out;
  }

  void Apply(std::function<void(TH &)> f) {
    f(CCInc);
    f(CC0Pi);
    f(CC1CPi);
    f(CC1Pi0);
    f(CCOther);
    f(NCInc);
    f(NC0Pi);
    f(NC1CPi);
    f(NC1Pi0);
    f(NCOther);
  }
};

TH1D *totxsecs;

hblob<TH1D> *Enu;
hblob<TH3D> *EnuPLepThetaLep;
hblob<TH3D> *EnuPLepEAvHad;
hblob<TH3D> *EnuQ2EAvHad;
hblob<TH3D> *EnuPtLepEAvHad;

void Fill(TTreeReader &rdr, bool ist2k, int tgta_select = 0) {
  Enu = new hblob<TH1D>(
      "Enu",
      ";#it{E_{#nu}} (GeV);#it{y}; d#sigma/d#it{E_{#nu}} cm^{2} GeV^{-1}", 200,
      0, 10);

  if (ist2k) {
    EnuPLepThetaLep = new hblob<TH3D>(
        "EnuPLepThetaLep",
        ";#it{E_{#nu}} (GeV);#it{p}_{lep} (GeV/#it{c});#theta_{lep} (rad)"
        "GeV^{-2} #it{c} rad^{-1}",
        100, 0, 5, 50, 0, 5, 25, 0, M_PI);
  } else {
    EnuPLepEAvHad = new hblob<TH3D>(
        "EnuPLepEAvHad",
        ";#it{E_{#nu}} (GeV);#it{p}_{lep} (GeV/#it{c}); #it{E}_{Had}^{Vis} "
        "(GeV)",
        100, 0, 5, 50, 0, 5, 40, 0, 4);
    EnuQ2EAvHad = new hblob<TH3D>(
        "EnuQ2EAvHad",
        ";#it{E_{#nu}} (GeV);Q^{2} (GeV^{2}/#it{c}^{2}); #it{E}_{Had}^{Vis} "
        "(GeV)",
        100, 0, 5, 50, 0, 5, 40, 0, 4);
    EnuPtLepEAvHad = new hblob<TH3D>(
        "EnuPtLepEAvHad",
        ";#it{E_{#nu}} (GeV);#it{p}^{T}_{lep} (GeV/#it{c}); #it{E}_{Had}^{Vis} "
        "(GeV)",
        100, 0, 5, 40, 0, 2, 40, 0, 4);
  }

  totxsecs = new TH1D("totxsecs", ";;#sigma^{#int#Phi} cm^{2}", 10, 0, 10);
  totxsecs->SetDirectory(nullptr);

  TTreeReaderValue<float> CosLep(rdr, "CosLep");
  TTreeReaderValue<float> EavAlt(rdr, "EavAlt");
  TTreeReaderValue<float> Enu_true(rdr, "Enu_true");
  TTreeReaderValue<float> PLep_v(rdr, "PLep");
  TTreeReaderValue<float> Q2(rdr, "Q2");
  TTreeReaderValue<int> tgta(rdr, "tgta");

  TTreeReaderValue<bool> flagCCINC(rdr, "flagCCINC");
  TTreeReaderValue<bool> flagCC0pi(rdr, "flagCC0pi");
  TTreeReaderValue<bool> flagCC1pip(rdr, "flagCC1pip");
  TTreeReaderValue<bool> flagCC1pim(rdr, "flagCC1pim");
  TTreeReaderValue<bool> flagCC1pi0(rdr, "flagCC1pi0");

  TTreeReaderValue<bool> flagNCINC(rdr, "flagNCINC");
  TTreeReaderValue<bool> flagNC0pi(rdr, "flagNC0pi");
  TTreeReaderValue<bool> flagNC1pip(rdr, "flagNC1pip");
  TTreeReaderValue<bool> flagNC1pim(rdr, "flagNC1pim");
  TTreeReaderValue<bool> flagNC1pi0(rdr, "flagNC1pi0");

  TTreeReaderValue<double> fScaleFactor(rdr, "fScaleFactor");
  TTreeReaderValue<double> RWWeight(rdr, "RWWeight");

  size_t nents = rdr.GetEntries(true);
  size_t ent_it = 0;
  size_t shout_it = nents / 25;

  while (rdr.Next()) {

    if (ent_it && !(ent_it % shout_it)) {
      std::cout << "[Read] " << ent_it << "/" << nents << "("
                << (100 * ent_it / nents) << "%)" << std::endl;
    }

    if (tgta_select && (tgta != tgta_select)) {
      ent_it++;
      continue;
    }

    double w = *fScaleFactor * *RWWeight;

    FlagBlob fblob{
        *flagCCINC, *flagCC0pi, *flagCC1pip || *flagCC1pim, *flagCC1pi0,
        *flagNCINC, *flagNC0pi, *flagNC1pip || *flagNC1pim, *flagNC1pi0};

    if (fblob.flagCCINC) {
      totxsecs->Fill(0.0, w);
      if (fblob.flagCC0pi) {
        totxsecs->Fill(1.0, w);
      } else if (fblob.flagCC1cpi) {
        totxsecs->Fill(2.0, w);
      } else if (fblob.flagCC1pi0) {
        totxsecs->Fill(3.0, w);
      } else {
        totxsecs->Fill(4.0, w);
      }
    } else if (fblob.flagNCINC) {
      totxsecs->Fill(5.0, w);
      if (fblob.flagNC0pi) {
        totxsecs->Fill(6.0, w);
      } else if (fblob.flagNC1cpi) {
        totxsecs->Fill(7.0, w);
      } else if (fblob.flagNC1pi0) {
        totxsecs->Fill(8.0, w);
      } else {
        totxsecs->Fill(9.0, w);
      }
    }

    Enu->Fill(w, fblob, *Enu_true);
    if (ist2k) {
      EnuPLepThetaLep->Fill(w, fblob, *Enu_true, *PLep_v, acos(*CosLep));
    } else {
      EnuPLepEAvHad->Fill(w, fblob, *Enu_true, *PLep_v, *EavAlt);
      EnuQ2EAvHad->Fill(w, fblob, *Enu_true, *Q2, *EavAlt);
      EnuPtLepEAvHad->Fill(w, fblob, *Enu_true,
                           (*PLep_v) * sqrt(1 - pow(*CosLep, 2)), *EavAlt);
    }
    ent_it++;
  }
}

int main(int argc, char const *argv[]) {
  if (argc < 4) {
    std::cout << "Expects 3 arguments." << std::endl;
    return 1;
  }

  gStyle->SetOptStat(false);

  TFile fin(argv[1]);
  if (fin.IsZombie()) {
    std::cout << "Failed to read " << argv[1] << std::endl;
    return 2;
  }

  TTreeReader rdr("T2KNOvATruthTree", &fin);

  bool ist2k = false;

  if (std::string(argv[2]) == "NOvAND") {
    ist2k = false;
  } else if (std::string(argv[2]) == "ND280") {
    ist2k = true;
  } else {
    std::cout << "[ERROR]: Expected option 2 to specify either NOvAND or ND280."
              << std::endl;
    return 1;
  }

  int tgta_select = 0;
  if (argc > 4) {
    if (std::string(argv[4]) == "C") {
      tgta_select = 12;
    } else if (std::string(argv[4]) == "H") {
      tgta_select = 1;
    } else if (std::string(argv[4]) == "O") {
      tgta_select = 6;
    } else if (std::string(argv[4]) == "any") {
      tgta_select = 0;
    } else {
      std::cout << "Invalid target selector passed: " <
          argv[4] << ". Should be C/H/O/any." return 1;
    }
  }

  Fill(rdr, ist2k, tgta_select);

  TFile fout(argv[3], "UPDATE");
  if (fin.IsZombie()) {
    std::cout << "Failed to write " << argv[3] << std::endl;
    return 2;
  }

  TDirectory *dout = &fout;
  if (argc > 5) {
    std::string fqdir = argv[5];
    while (fqdir.find("/") != std::string::npos) {
      std::string ndir = fqdir.substr(0, fqdir.find("/"));
      if (!dout->GetDirectory(ndir.c_str())) {
        dout = dout->mkdir(ndir.c_str());
      } else {
        dout = dout->GetDirectory(ndir.c_str());
      }
      fqdir = fqdir.substr(fqdir.find("/") + 1);
    }
    dout = dout->mkdir(fqdir.c_str());
  }

  totxsecs->GetXaxis()->SetBinLabel(1, "CCInc");
  totxsecs->GetXaxis()->SetBinLabel(2, "CC0#pi");
  totxsecs->GetXaxis()->SetBinLabel(3, "CC1#pi^{#pm}");
  totxsecs->GetXaxis()->SetBinLabel(4, "CC1#pi^{0}");
  totxsecs->GetXaxis()->SetBinLabel(5, "CCMulti#pi + CCOther");
  totxsecs->GetXaxis()->SetBinLabel(6, "NCInc");
  totxsecs->GetXaxis()->SetBinLabel(7, "NC0#pi");
  totxsecs->GetXaxis()->SetBinLabel(8, "NC1#pi^{#pm}");
  totxsecs->GetXaxis()->SetBinLabel(9, "NC1#pi^{0}");
  totxsecs->GetXaxis()->SetBinLabel(10, "NCMulti#pi + CCOther");
  dout->WriteTObject(totxsecs, "totxsecs");

  Enu->Write(dout);

  auto ProjYZ = [=](TH3D const &h) -> TH2D {
    std::unique_ptr<TH2D> h2(dynamic_cast<TH2D *>(h.Project3D("yz")));
    h2->SetDirectory(nullptr);
    h2->SetName((std::string(h.GetName()) + "_pyz").c_str());
    return TH2D(*h2.get());
  };
  auto ProjYX = [=](TH3D const &h) -> TH2D {
    std::unique_ptr<TH2D> h2(dynamic_cast<TH2D *>(h.Project3D("yx")));
    h2->SetDirectory(nullptr);
    h2->SetName((std::string(h.GetName()) + "_pyx").c_str());
    return TH2D(*h2.get());
  };
  auto ProjXZ = [=](TH3D const &h) -> TH2D {
    std::unique_ptr<TH2D> h2(dynamic_cast<TH2D *>(h.Project3D("xz")));
    h2->SetDirectory(nullptr);
    h2->SetName((std::string(h.GetName()) + "_pxz").c_str());
    return TH2D(*h2.get());
  };

  if (ist2k) {
    EnuPLepThetaLep->Write(dout);
    auto PLepThetaLep = EnuPLepThetaLep->Transform<TH2D>(ProjYZ);
    PLepThetaLep.Write(dout);
    auto EnuPLep = EnuPLepThetaLep->Transform<TH2D>(ProjYX);
    EnuPLep.Write(dout);
    auto EnuThetaLep = EnuPLepThetaLep->Transform<TH2D>(ProjXZ);
    EnuThetaLep.Write(dout);

  } else {
    EnuPLepEAvHad->Write(dout);
    auto PLepEAvHad = EnuPLepEAvHad->Transform<TH2D>(ProjYZ);
    PLepEAvHad.Write(dout);
    auto EnuPLep = EnuPLepEAvHad->Transform<TH2D>(ProjYX);
    EnuPLep.Write(dout);
    auto EnuEAvHad = EnuPLepEAvHad->Transform<TH2D>(ProjXZ);
    EnuEAvHad.Write(dout);

    EnuQ2EAvHad->Write(dout);
    auto Q2EAvHad = EnuQ2EAvHad->Transform<TH2D>(ProjYZ);
    Q2EAvHad.Write(dout);
    auto EnuQ2 = EnuQ2EAvHad->Transform<TH2D>(ProjYX);
    EnuQ2.Write(dout);

    EnuPtLepEAvHad->Write(dout);
    auto PtLepEAvHad = EnuPtLepEAvHad->Transform<TH2D>(ProjYZ);
    PtLepEAvHad.Write(dout);
    auto EnuPtLep = EnuPtLepEAvHad->Transform<TH2D>(ProjYX);
    EnuPtLep.Write(dout);
  }

  fout.Close();
}
