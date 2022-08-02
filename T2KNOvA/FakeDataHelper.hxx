#pragma once

#include <cmath>
#include <iostream>
#include <sstream>
#include <unordered_map>

#include "T2KNOvA/ROOTHelper.hxx"
#include "T2KNOvA/TrueSelectionHelper.hxx"

#include "TDirectory.h"
#include "TFile.h"
#include "TH3.h"
#include "THn.h"

// #define DEBUG_HIST_LOADER

namespace t2knova {

double GetFakeDataWeight_NOvAToT2KND_PtLep(int nu_pdg, int lep_pdg, int tgta,
                                           double E_nu_GeV, double PtLep_GeV,
                                           double EVisHadronic_GeV,
                                           int PrimSel = -1);
double GetFakeDataWeight_NOvAToT2KPre_PtLep(int nu_pdg, int lep_pdg, int tgta,
                                            double E_nu_GeV, double PtLep_GeV,
                                            double EVisHadronic_GeV,
                                            int PrimSel = -1);
double GetFakeDataWeight_NOvAToT2KMnv1Pi_PtLep(int nu_pdg, int lep_pdg,
                                               int tgta, double E_nu_GeV,
                                               double PtLep_GeV,
                                               double EVisHadronic_GeV,
                                               int PrimSel = -1);
double GetFakeDataWeight_NOvAToT2KNonQE_PtLep(int nu_pdg, int lep_pdg, int tgta,
                                              double E_nu_GeV, double PtLep_GeV,
                                              double EVisHadronic_GeV,
                                              int PrimSel = -1);
double GetFakeDataWeight_ND280ToNOvA(int nu_pdg, int lep_pdg, int tgta,
                                     double E_nu_GeV, double PLep_GeV,
                                     double ThetaLep, int PrimSel = -1);

double GetFakeDataWeight_ND280ToT2KNonQE(int nu_pdg, int lep_pdg, int tgta,
                                         double E_nu_GeV, double PLep_GeV,
                                         double ThetaLep, int PrimSel = -1);
} // namespace t2knova

namespace t2knova {

enum reweightconfig {
  kNoWeight = 0,
  kT2KND_to_NOvA,
  kT2KND_to_T2KNonQE,
  kNOvA_to_T2KND_ptlep,
  kNOvA_to_T2KPre_ptlep,
  kNOvA_to_T2KMnv1Pi_ptlep,
  kNOvA_to_T2KNonQE_ptlep,
};

const char *all_tgta_str[] = {"H", "C", "O"};
const int all_tgta[] = {1, 12, 16};

static bool loaded = false;
static std::unordered_map<
    int, std::unordered_map<int, std::unordered_map<int, std::unique_ptr<TH1>>>>
    rwhists;

template <typename T, size_t N> inline size_t arrsize(T (&arr)[N]) { return N; }
template <typename T> inline size_t arrsize(T arr) { return arr.size(); }

inline void LoadHists(std::string const &inputfile,
                      std::map<reweightconfig, std::string> const &inputhists) {
  std::unique_ptr<TFile> fin(new TFile(inputfile.c_str()));
  if (fin->IsZombie() || !fin->IsOpen()) {
    std::cout << "Failed to read \"" << inputfile << "\"" << std::endl;
    abort();
  }
  int found = 0;
  for (size_t i = 0; i < arrsize(all_nuspecies); ++i) {
    nuspecies nuspec = nuspecies(i);
    std::string nuspec_str = all_nuspecies[i];

    for (auto &rwconfighistname : inputhists) {
      reweightconfig rwconfig = rwconfighistname.first;

      std::string rwconfig_str = rwconfighistname.second;

      for (selection sel : ReWeightSelectionList) {
        std::string sel_str = SelectionList[sel];

        for (size_t l = 0; l < arrsize(all_tgta_str); ++l) {
          std::string tgta_str = all_tgta_str[l];
          int tgta_sel_offset = all_tgta[l] * 100;

          rwhists[nuspec][rwconfig][tgta_sel_offset + sel] =
              GetTH1(fin, rwconfig_str + "/" + tgta_str + "/" + nuspec_str +
                              "/" + sel_str);
          if (rwhists[nuspec][rwconfig][tgta_sel_offset + sel]) {
            rwhists[nuspec][rwconfig][tgta_sel_offset + sel]->SetDirectory(
                nullptr);
            found++;
#ifdef DEBUG_HIST_LOADER
            std::cout << "[INFO]: Loaded: "
                      << rwconfig_str + "/" + tgta_str + "/" + nuspec_str +
                             "/" + sel_str
                      << std::endl;
#endif
          }
#ifdef DEBUG_HIST_LOADER
          else {
            std::cout << "[WARN]: Failed to load: "
                      << rwconfig_str + "/" + tgta_str + "/" + nuspec_str +
                             "/" + sel_str
                      << std::endl;
          }
#endif
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

inline double
GetFakeDataWeight_NOvAToT2K_PtLep(reweightconfig conf, int nu_pdg, int lep_pdg,
                                  int tgta, double E_nu_GeV, double PtLep_GeV,
                                  double EVisHadronic_GeV, int PrimSel) {
  if (!loaded) {
    std::cout << "[ERROR]: Have not loaded t2knova reweight histogram ratios."
              << std::endl;
    abort();
  }

  if (PrimSel == kNoPrimarySel) {
    return 1;
  }

  nuspecies nuspec = getnuspec(nu_pdg);

  bool iscc = (nu_pdg != lep_pdg);

  int offset = (PrimSel == -1) ? (iscc ? kCCInc : kNCInc) : PrimSel;

  std::unique_ptr<TH1> &rathist = rwhists[nuspec][conf][(tgta * 100) + offset];

  if (rathist) {
    double wght =
        EvalHist3D(rathist, E_nu_GeV, PtLep_GeV, EVisHadronic_GeV, false);
    return wght && !std::isnormal(wght) ? 0 : wght;
  }
  return 1;
}

inline double GetFakeDataWeight_NOvAToT2KND_PtLep(int nu_pdg, int lep_pdg,
                                                  int tgta, double E_nu_GeV,
                                                  double PtLep_GeV,
                                                  double EVisHadronic_GeV,
                                                  int PrimSel) {
  return GetFakeDataWeight_NOvAToT2K_PtLep(kNOvA_to_T2KND_ptlep, nu_pdg,
                                           lep_pdg, tgta, E_nu_GeV, PtLep_GeV,
                                           EVisHadronic_GeV, PrimSel);
}

inline double GetFakeDataWeight_NOvAToT2KPre_PtLep(int nu_pdg, int lep_pdg,
                                                   int tgta, double E_nu_GeV,
                                                   double PtLep_GeV,
                                                   double EVisHadronic_GeV,
                                                   int PrimSel) {
  return GetFakeDataWeight_NOvAToT2K_PtLep(kNOvA_to_T2KPre_ptlep, nu_pdg,
                                           lep_pdg, tgta, E_nu_GeV, PtLep_GeV,
                                           EVisHadronic_GeV, PrimSel);
}

inline double GetFakeDataWeight_NOvAToT2KMnv1Pi_PtLep(int nu_pdg, int lep_pdg,
                                                      int tgta, double E_nu_GeV,
                                                      double PtLep_GeV,
                                                      double EVisHadronic_GeV,
                                                      int PrimSel) {
  return GetFakeDataWeight_NOvAToT2K_PtLep(kNOvA_to_T2KMnv1Pi_ptlep, nu_pdg,
                                           lep_pdg, tgta, E_nu_GeV, PtLep_GeV,
                                           EVisHadronic_GeV, PrimSel);
}

inline double GetFakeDataWeight_NOvAToT2KNonQE_PtLep(int nu_pdg, int lep_pdg,
                                                     int tgta, double E_nu_GeV,
                                                     double PtLep_GeV,
                                                     double EVisHadronic_GeV,
                                                     int PrimSel) {
  return GetFakeDataWeight_NOvAToT2K_PtLep(kNOvA_to_T2KNonQE_ptlep, nu_pdg,
                                           lep_pdg, tgta, E_nu_GeV, PtLep_GeV,
                                           EVisHadronic_GeV, PrimSel);
}

inline double GetFakeDataWeight_ND280ToNOvA(int nu_pdg, int lep_pdg, int tgta,
                                            double E_nu_GeV, double PLep_GeV,
                                            double ThetaLep, int PrimSel) {
  if (!loaded) {
    std::cout << "[ERROR]: Have not loaded t2knova reweight histogram ratios."
              << std::endl;
    abort();
  }

  if (PrimSel == kNoPrimarySel) {
    return 1;
  }

  nuspecies nuspec = getnuspec(nu_pdg);

  std::unique_ptr<TH1> &rathist =
      rwhists[nuspec][kT2KND_to_NOvA][(tgta * 100) + PrimSel];

  if (rathist) {
    return EvalHist3D(rathist, E_nu_GeV, PLep_GeV, ThetaLep, false);
  }

  return 1;
}

inline double GetFakeDataWeight_ND280ToT2KNonQE(int nu_pdg, int lep_pdg,
                                                int tgta, double E_nu_GeV,
                                                double PLep_GeV,
                                                double ThetaLep, int PrimSel) {
  if (!loaded) {
    std::cout << "[ERROR]: Have not loaded t2knova reweight histogram ratios."
              << std::endl;
    abort();
  }

  if (PrimSel == kNoPrimarySel) {
    return 1;
  }

  nuspecies nuspec = getnuspec(nu_pdg);

  std::unique_ptr<TH1> &rathist =
      rwhists[nuspec][kT2KND_to_T2KNonQE][(tgta * 100) + PrimSel];

  if (rathist) {
    return EvalHist3D(rathist, E_nu_GeV, PLep_GeV, ThetaLep, false);
  }

  return 1;
}


///\note See https://arxiv.org/abs/1903.01558
inline double GetMINERvASPPLowQ2SuppressionWeight(double Q2_True_GeV,
                                                  double parameter_value = 1) {

  // Fit parameters for FrInel+Low-Q2 tune were
  // MA_RES = 0.93 +/- 0.05
  // NormRes = 116 +/- 7
  // NonRes1pi = 46 +/- 4
  // NonRes2pi = 120 +/- 32
  // ThetaPi = 1.0 (at limit)
  // FrInel = 132 +/- 27
  // R1 = 0.37 +/- 0.09
  // R2 = 0.60 +/- 0.16
  //
  // Fit parameters for FrAbs+Low-Q2 tune
  // MA_RES = 0.92 +/- 0.02
  // NormRes = 116 +/- 3
  // NonRes1pi = 46 +/- 4
  // NonRes2pi = 99 +/- 31
  // ThetaPi = 1.0 (at limit)
  // FrAbs = 48 +/- 21
  // R1 = 0.32 +/- 0.06
  // R2 = 0.5 (limit)
  static double const Q2_Max = 0.7;
  static double const Q2_t1 = 0;
  static double const Q2_t2 = 0.35;
  static double const R1 = 0.37;
  static double const R2 = 0.6;

  if ((Q2_True_GeV > Q2_Max) || (Q2_True_GeV < 0)) {
    return 1;
  }

  double RQ2 = (R2 * ((Q2_True_GeV - Q2_t1) * (Q2_True_GeV - Q2_Max)) /
                ((Q2_t2 - Q2_t1) * (Q2_t2 - Q2_Max))) +
               (((Q2_True_GeV - Q2_t1) * (Q2_True_GeV - Q2_t2)) /
                ((Q2_Max - Q2_t1) * (Q2_Max - Q2_t2)));
  return 1 - parameter_value * ((1 - R1) * pow((1 - RQ2), 2));
}

inline double UnWeightQ2BinWeights_T2K2020(double Q2_True_GeV) {
  static double Q2Weights[] = {0.7841, 0.8868, 1.0228, 1.0268,
                               1.0867, 1.2568, 1.1360, 1.2593};
  static double fBinWidth_GeV2 = 0.05;
  double Q2Weight = 1;
  if (Q2_True_GeV <= 0.25) {
    Q2Weight = Q2Weights[int(std::floor(Q2_True_GeV / fBinWidth_GeV2))];
  } else if (Q2_True_GeV > 0.25 && Q2_True_GeV <= 0.5) {
    Q2Weight = Q2Weights[5];
  } else if (Q2_True_GeV > 0.5 && Q2_True_GeV <= 1.0) {
    Q2Weight = Q2Weights[6];
  } else if (Q2_True_GeV > 1.0) {
    Q2Weight = Q2Weights[7];
  }
  return 1.0 / Q2Weight;
}

inline double GetnonQEWeight(int nuPDG, double Q2_Reco_GeV) {

  static bool first = true;
  static TH1 *nuWeights = nullptr;
  static TH1 *nubWeights = nullptr;

  if (first) {
    char *T2KNOVA_INPUTS = getenv("T2KNOVA_INPUTS");

    if (!T2KNOVA_INPUTS) {
      std::cout << "[ERROR]: Expected T2KNOVA_INPUTS environment variable to "
                   "be defined."
                << std::endl;
      abort();
    }

    TFile *fin = new TFile(
        (std::string(T2KNOVA_INPUTS) + "/ScalingHisto_nu_antinu.root").c_str());
    if (!fin || !fin->IsOpen()) {
      std::cout << "[ERROR]; Failed to open " << T2KNOVA_INPUTS
                << "/ScalingHisto_nu_antinu.root" << std::endl;
      abort();
    }

    fin->GetObject("ScalingHisto_FGD1_numu", nuWeights);
    fin->GetObject("ScalingHisto_FGD1_anumu", nubWeights);

    nuWeights = dynamic_cast<TH1 *>(nuWeights->Clone("ScalingHisto_FGD1_numu"));
    if (!nubWeights) {
      std::cout << "[ERROR]: Failed to read ScalingHisto_FGD1_numu from "
                << T2KNOVA_INPUTS << "/ScalingHisto_nu_antinu.root"
                << std::endl;
      abort();
    }
    nuWeights->SetDirectory(nullptr);
    nubWeights =
        dynamic_cast<TH1 *>(nubWeights->Clone("ScalingHisto_FGD1_anumu"));
    if (!nubWeights) {
      std::cout << "[ERROR]: Failed to read ScalingHisto_FGD1_anumu from "
                << T2KNOVA_INPUTS << "/ScalingHisto_nu_antinu.root"
                << std::endl;
      abort();
    }
    nubWeights->SetDirectory(nullptr);

    fin->Close();
    delete fin;
    first = false;
  }

  TH1 *h = nuPDG > 0 ? nuWeights : nubWeights;

  int bin = h->FindFixBin(Q2_Reco_GeV);

  return h->GetBinContent(bin);
}

} // namespace t2knova
