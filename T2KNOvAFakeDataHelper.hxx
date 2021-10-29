#pragma once

#include "TDirectory.h"
#include "TFile.h"
#include "TH3.h"
#include "THn.h"

#include <iostream>
#include <unordered_map>

namespace t2knova {

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

struct FSParticleSummary {
  FSParticleSummary()
      : IsCC(false), NChLep(0), NNeutralLep(0), NProton(0), NNeutron(0),
        NAntiNucleon(0), NChPi(0), NPi0(0), NGamma(0), NKaons(0), NOther(0),
        NNuclear(0) {}

  bool IsCC;
  int NChLep;
  int NNeutralLep;
  int NProton;
  int NNeutron;
  int NAntiNucleon;
  int NChPi;
  int NPi0;
  int NGamma;
  int NKaons;
  int NOther;
  int NNuclear;

  int GetNNonNucleonPions() {
    int count = 0;

    int NChLep_extra = NChLep - (IsCC ? 1 : 0);
    int NNeutralLep_extra = NNeutralLep - (IsCC ? 0 : 1);

    if ((NChLep_extra < 0) || (NNeutralLep_extra < 0)) {
      std::cout << "Found bad event: NChLep_extra = " << NChLep_extra
                << ", NNeutralLep_extra = " << NNeutralLep_extra << std::endl;
      abort();
    }

    count += NChLep_extra;
    count += NNeutralLep_extra;
    count += NAntiNucleon;
    count += NGamma;
    count += NKaons;
    count += NOther;

    return count;
  }
};

inline FSParticleSummary T2KNOvAFlatTreeToFSParticleSummary(int NFSP,
                                                            int *FSPDG) {

  int PDG_FSLep = 0;

  FSParticleSummary fsps;

  for (int i = 0; i < NFSP; ++i) {
    switch (FSPDG[i]) {
    case 11:
    case 13:
    case 15:
    case -11:
    case -13:
    case -15: {
      if (!PDG_FSLep) {
        PDG_FSLep = FSPDG[i];
      }
      fsps.NChLep++;
      break;
    }

    case 12:
    case 14:
    case 16:
    case -12:
    case -14:
    case -16: {
      if (!PDG_FSLep) {
        PDG_FSLep = FSPDG[i];
      }
      fsps.NNeutralLep++;

      break;
    }
    case 2212: {
      fsps.NProton++;
      break;
    }
    case 2112: {
      fsps.NNeutron++;
      break;
    }
    case -2212:
    case -2112: {
      fsps.NAntiNucleon++;
      break;
    }
    case 211:
    case -211: {
      fsps.NChPi++;
      break;
    }
    case 111: {
      fsps.NPi0++;
      break;
    }
    case 22: {
      fsps.NGamma++;
      break;
    }
    case 321:
    case 311:
    case -321:
    case -311: {
      fsps.NKaons++;
      break;
    }

    default: {
      if (FSPDG[i] > 1000000000) {
        fsps.NNuclear++
      } else {
        fsps.NOther++
      }
      break;
    }
    }
  }
  fsps.IsCC = PDG_FSLep % 2;

  return fsps;
}

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
