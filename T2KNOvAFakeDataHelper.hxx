#pragma once

#include "TDirectory.h"
#include "TFile.h"
#include "TH3.h"
#include "THn.h"

#include <iostream>
#include <unordered_map>

namespace t2knova {

struct FlagBlob {
  FlagBlob()
      : Mode(0), flagCCINC(false), flagCC0pi(false), flagCC1cpi(false),
        flagCC1pi0(false), flagNCINC(false), flagNC0pi(false),
        flagNC1cpi(false), flagNC1pi0(false) {}
  int Mode;
  bool flagCCINC;
  bool flagCC0pi;
  bool flagCC1cpi;
  bool flagCC1pi0;
  bool flagNCINC;
  bool flagNC0pi;
  bool flagNC1cpi;
  bool flagNC1pi0;
};

inline FlagBlob GetFlagBlob(bool iscc, int NFSCpi, int NFSpi0, int NFSOther,
                            int Mode = 0) {

  FlagBlob fb;
  fb.Mode = Mode;

  fb.flagCCINC = iscc;
  fb.flagNCINC = !iscc;

  if (NFSOther == 0) { // Normalish event with just nucleons and pions
    if ((NFSpi0 + NFSCPi) == 0) {
      fb.flagCC0pi = iscc;
      fb.flagNC0pi = !iscc;
    } else if ((NFSpi0 == 0) && ((NFSCPi) == 1)) {
      fb.flagCC1cpi = iscc;
      fb.flagNC1cpi = !iscc;
    } else if ((NFSpi0 == 1) && ((NFSCPi) == 0)) {
      fb.flagCC1pi0 = iscc;
      fb.flagNC1pi0 = !iscc;
    } // else CCOther
  }   // else CCOther

  return fb;
}

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
  std::unordered_map<int, TH> ModeHists;

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

    for (auto &h : ModeHists) {
      h.second.SetName((name + "_Mode_" + (h.first < 0 ? "m" : "") +
                        std::to_string(std::abs(h.first)))
                           .c_str());
    }
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

    if (blb.Mode != 0) {
      if (!ModeHists.count(blb.Mode)) {
        ModeHists.emplace(blb.Mode, CCInc);

        ModeHists[blb.Mode].Reset();

        ModeHists[blb.Mode].SetDirectory(nullptr);
      }
      ModeHists[blb.Mode].Fill(xy..., w);
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

    for (auto const &h : ModeHists) {
      out.ModeHists[h.first] = THOut(f(h.second));
    }
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

    for (auto &h : ModeHists) {
      f(h.second);
    }
  }
};

double GetFakeDataWeight_NOvAToT2K_PLep(int nu_pdg, int lep_pdg, int tgta,
                                        double E_nu_GeV, double PLep_GeV,
                                        double EVisHadronic_GeV,
                                        bool interpolate = true);
double GetFakeDataWeight_NOvAToT2K_Q2(int nu_pdg, int lep_pdg, int tgta,
                                      double E_nu_GeV, double Q2_GeV2,
                                      double EVisHadronic_GeV,
                                      bool interpolate = true);
double GetFakeDataWeight_NOvAToT2K_PtLep(int nu_pdg, int lep_pdg, int tgta,
                                         double E_nu_GeV, double PtLep_GeV,
                                         double EVisHadronic_GeV,
                                         bool interpolate = true);
double GetFakeDataWeight_ND280ToNOvA(int nu_pdg, int lep_pdg, int tgta,
                                     double E_nu_GeV, double PLep_GeV,
                                     double ThetaLep, int NFSCPi, int NFSPi0,
                                     int NOther, bool interpolate = true);

double GetFakeDataWeight_ND280ToNOvA_Enu(int nu_pdg, int lep_pdg, int tgta,
                                         double E_nu_GeV, int NFSCPi,
                                         int NFSPi0, int NOther,
                                         bool interpolate = true);

double GetFakeDataWeight_ND280ToNOvA_Q2(int nu_pdg, int lep_pdg, int tgta,
                                        double Q2_GeV2, int NFSCPi, int NFSPi0,
                                        int NOther, bool interpolate = true);
} // namespace t2knova

namespace t2knova {

inline TH1 *GetTH1(TFile *f, std::string const &name) {
  TDirectory *odir = gDirectory;

  TH1 *h;
  f->GetObject(name.c_str(), h);
  if (h) {
    h->SetDirectory(nullptr);
  }

  if (odir) {
    gDirectory->cd();
  }

  return h;
} // namespace t2knova

enum nuspecies { kNuMu = 0, kNuMub, kNuE, kNuEb };
const char *all_nuspecies[] = {"numu", "numub", "nue", "nueb"};
const char *all_nuspecies_latex[] = {"#nu_{#mu}", "#bar{#nu}_{#mu}", "#nu_{e}",
                                     "#bar{#nu}_{e}"};

inline nuspecies getnuspec(int pdg) {
  switch (pdg) {
  case 14: {
    return kNuMu;
  }
  case -14: {
    return kNuMub;
  }
  case 12: {
    return kNuE;
  }
  case -12: {
    return kNuEb;
  }
  default: {
    std::cerr << "[ERROR]: Invalid neutrino PDG " << pdg
              << " passed to t2knova::getnuspec" << std::endl;
    abort();
  }
  }
}

enum reweightconfig {
  kT2KND_to_NOvA = 0,
  kT2KND_to_NOvA_Enu,
  kT2KND_to_NOvA_Q2,
  kNOvA_to_T2KND_plep,
  kNOvA_to_T2KND_Q2,
  kNOvA_to_T2KND_ptlep,
  kNoWeight
};
const char *all_rwconfig[] = {"t2knd_to_nova",    "t2knd_to_nova_Enu",
                              "t2knd_to_nova_Q2", "nova_to_t2k_plep",
                              "nova_to_t2k_Q2",   "nova_to_t2k_ptlep"};

enum selection {
  kCCINC = 0,
  kCC0pi,
  kCC1cpi,
  kCC1pi0,
  kCCOther,
  kNCINC,
  kNC0pi,
  kNC1cpi,
  kNC1pi0,
  kNCOther
};
const char *all_sel[] = {"CCInc", "CC0Pi", "CC1CPi", "CC1Pi0", "CCOther",
                         "NCInc", "NC0Pi", "NC1CPi", "NC1Pi0", "NCOther"};

const char *all_tgta_str[] = {"H", "C", "O"};
const int all_tgta[] = {1, 12, 16};

static bool loaded = false;
static std::unordered_map<
    nuspecies,
    std::unordered_map<reweightconfig, std::unordered_map<int, TH1 *>>>
    rwhists;

template <typename T, size_t N> inline size_t arrsize(T (&arr)[N]) { return N; }

inline void LoadHists(std::string const &inputfile = "FakeDataInputs.root") {
  TFile fin(inputfile.c_str());
  if (fin.IsZombie()) {
    std::cout << "Failed to read \"" << inputfile << "\"" << std::endl;
    abort();
  }
  int found = 0;
  for (int i = 0; i < arrsize(all_nuspecies); ++i) {
    nuspecies nuspec = nuspecies(i);
    std::string nuspec_str = all_nuspecies[i];

    for (int j = 0; j < arrsize(all_rwconfig); ++j) {
      reweightconfig rwconfig = reweightconfig(j);
      std::string rwconfig_str = all_rwconfig[j];

      for (int k = 0; k < arrsize(all_sel); ++k) {
        selection sel = selection(k);
        std::string sel_str = all_sel[k];

        for (int l = 0; l < arrsize(all_tgta_str); ++l) {
          std::string tgta_str = all_tgta_str[l];
          int tgta_sel_offset = all_tgta[l] * 100;

          rwhists[nuspec][rwconfig][tgta_sel_offset + sel] =
              dynamic_cast<TH1 *>(
                  t2knova::GetTH1(&fin, rwconfig_str + "_" + tgta_str + "_" +
                                            nuspec_str + "_" + sel_str));
          if (rwhists[nuspec][rwconfig][tgta_sel_offset + sel]) {
            rwhists[nuspec][rwconfig][tgta_sel_offset + sel]->SetDirectory(
                nullptr);
            found++;
          } else {
            std::cout << "[WARN]: Couldn't find expected histogram: "
                      << (rwconfig_str + "_" + tgta_str + "_" + nuspec_str +
                          "_" + sel_str)
                      << std::endl;
          }
        }
      }
    }
  }
  if (found == 0) {
    std::cout << "[ERROR]: Failed to find any histograms." << std::endl;
    abort();
  }
  loaded = true;
}

inline double EvalHist3D(TH1 *h1, double x, double y, double z,
                         bool interpolate = true) {
  TH3 *h = static_cast<TH3 *>(h1);

  int xbin = h->GetXaxis()->FindFixBin(x);
  int ybin = h->GetYaxis()->FindFixBin(y);
  int zbin = h->GetZaxis()->FindFixBin(z);

  if ((xbin != 0) && (xbin != (h->GetXaxis()->GetNbins() + 1)) && (ybin != 0) &&
      (ybin != (h->GetYaxis()->GetNbins() + 1)) && (zbin != 0) &&
      (zbin != (h->GetZaxis()->GetNbins() + 1))) {

    if (interpolate && (xbin > 1) && (xbin < h->GetXaxis()->GetNbins()) &&
        (ybin > 1) && (ybin < h->GetYaxis()->GetNbins()) && (zbin > 1) &&
        (zbin < h->GetZaxis()->GetNbins())) {
      return h->Interpolate(x, y, z);
    }
    return h->GetBinContent(xbin, ybin, zbin);
  }

  return 1;
}

inline double EvalHist1D(TH1 *h, double x, bool interpolate = true) {

  int xbin = h->GetXaxis()->FindFixBin(x);

  if ((xbin != 0) && (xbin != (h->GetXaxis()->GetNbins() + 1))) {

    if (interpolate && (xbin > 1) && (xbin < h->GetXaxis()->GetNbins())) {
      return h->Interpolate(x);
    }
    return h->GetBinContent(xbin);
  }

  return 1;
}

inline double GetFakeDataWeight_NOvAToT2K_PLep(int nu_pdg, int lep_pdg,
                                               int tgta, double E_nu_GeV,
                                               double PLep_GeV,
                                               double EVisHadronic_GeV,
                                               bool interpolate) {
  if (!loaded) {
    LoadHists();
  }
  nuspecies nuspec = getnuspec(nu_pdg);

  bool iscc = (nu_pdg != lep_pdg);

  TH1 *rathist = rwhists[nuspec][kNOvA_to_T2KND_plep]
                        [(tgta * 100) + (iscc ? kCCINC : kNCINC)];
  if (rathist) {
    return EvalHist3D(rathist, E_nu_GeV, PLep_GeV, EVisHadronic_GeV,
                      interpolate);
  }

  return 1;
}

inline double GetFakeDataWeight_NOvAToT2K_Q2(int nu_pdg, int lep_pdg, int tgta,
                                             double E_nu_GeV, double Q2_GeV2,
                                             double EVisHadronic_GeV,
                                             bool interpolate) {
  if (!loaded) {
    LoadHists();
  }
  nuspecies nuspec = getnuspec(nu_pdg);

  bool iscc = (nu_pdg != lep_pdg);

  TH1 *rathist = rwhists[nuspec][kNOvA_to_T2KND_Q2]
                        [(tgta * 100) + (iscc ? kCCINC : kNCINC)];
  if (rathist) {
    return EvalHist3D(rathist, E_nu_GeV, Q2_GeV2, EVisHadronic_GeV,
                      interpolate);
  }
  return 1;
}

inline double GetFakeDataWeight_NOvAToT2K_PtLep(int nu_pdg, int lep_pdg,
                                                int tgta, double E_nu_GeV,
                                                double PtLep_GeV,
                                                double EVisHadronic_GeV,
                                                bool interpolate) {
  if (!loaded) {
    LoadHists();
  }
  nuspecies nuspec = getnuspec(nu_pdg);

  bool iscc = (nu_pdg != lep_pdg);

  TH1 *rathist = rwhists[nuspec][kNOvA_to_T2KND_ptlep]
                        [(tgta * 100) + (iscc ? kCCINC : kNCINC)];
  if (rathist) {
    return EvalHist3D(rathist, E_nu_GeV, PtLep_GeV, EVisHadronic_GeV,
                      interpolate);
  }
  return 1;
}

selection gett2ksel(bool iscc, int NFSCPi, int NFSPi0, int NOther) {
  selection sel = kCCINC;
  if ((NFSCPi + NFSPi0 + NOther) == 0) {
    sel = iscc ? kCC0pi : kNC0pi;
  } else if ((NFSCPi + NFSPi0 + NOther) == 1) {
    if (NFSCPi == 1) {
      sel = iscc ? kCC1cpi : kNC1cpi;
    } else if (NFSPi0 == 1) {
      sel = iscc ? kCC1pi0 : kNC1pi0;
    } else { // probably a kaon
      sel = iscc ? kCCOther : kNCOther;
    }
  } else {
    sel = iscc ? kCCOther : kNCOther;
  }
  return sel;
}

inline double GetFakeDataWeight_ND280ToNOvA(int nu_pdg, int lep_pdg, int tgta,
                                            double E_nu_GeV, double PLep_GeV,
                                            double ThetaLep, int NFSCPi,
                                            int NFSPi0, int NOther,
                                            bool interpolate) {
  if (!loaded) {
    LoadHists();
  }
  nuspecies nuspec = getnuspec(nu_pdg);

  bool iscc = (nu_pdg != lep_pdg);

  selection sel = gett2ksel(iscc, NFSCPi, NFSPi0, NOther);

  TH1 *rathist = rwhists[nuspec][kT2KND_to_NOvA][(tgta * 100) + sel];
  if (rathist) {
    return EvalHist3D(rathist, E_nu_GeV, PLep_GeV, ThetaLep, interpolate);
  }

  return 1;
}

inline double GetFakeDataWeight_ND280ToNOvA_Enu(int nu_pdg, int lep_pdg,
                                                int tgta, double E_nu_GeV,
                                                int NFSCPi, int NFSPi0,
                                                int NOther, bool interpolate) {
  if (!loaded) {
    LoadHists();
  }
  nuspecies nuspec = getnuspec(nu_pdg);

  bool iscc = (nu_pdg != lep_pdg);

  selection sel = gett2ksel(iscc, NFSCPi, NFSPi0, NOther);

  TH1 *rathist = rwhists[nuspec][kT2KND_to_NOvA_Enu][(tgta * 100) + sel];
  if (rathist) {
    return EvalHist1D(rathist, E_nu_GeV, interpolate);
  }

  return 1;
}

inline double GetFakeDataWeight_ND280ToNOvA_Q2(int nu_pdg, int lep_pdg,
                                               int tgta, double Q2_GeV2,
                                               int NFSCPi, int NFSPi0,
                                               int NOther, bool interpolate) {
  if (!loaded) {
    LoadHists();
  }
  nuspecies nuspec = getnuspec(nu_pdg);

  bool iscc = (nu_pdg != lep_pdg);

  selection sel = gett2ksel(iscc, NFSCPi, NFSPi0, NOther);

  TH1 *rathist = rwhists[nuspec][kT2KND_to_NOvA_Q2][(tgta * 100) + sel];
  if (rathist) {
    return EvalHist1D(rathist, Q2_GeV2, interpolate);
  }

  return 1;
}

} // namespace t2knova
