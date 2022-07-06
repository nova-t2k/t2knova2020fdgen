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
} // namespace t2knova

namespace t2knova {

enum reweightconfig {
  kT2KND_to_NOvA = 0,
  kNOvA_to_T2KND_ptlep,
  kNOvA_to_T2KPre_ptlep,
  kNOvA_to_T2KMnv1Pi_ptlep,
  kNOvA_to_T2KNonQE_ptlep,
  kNoWeight
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

} // namespace t2knova
