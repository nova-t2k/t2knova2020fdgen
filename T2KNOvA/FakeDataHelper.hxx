#pragma once

#include "T2KNOvA/ROOTHelper.hxx"
#include "T2KNOvA/TrueSelectionHelper.hxx"

#include "TDirectory.h"
#include "TFile.h"
#include "TH3.h"
#include "THn.h"

#include <iostream>
#include <sstream>
#include <unordered_map>

namespace t2knova {

constexpr double NMinEvs = 1;
constexpr double MaxFracError = 1.0 / sqrt(NMinEvs);

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
                                     double ThetaLep, int PrimSel,
                                     bool interpolate = true);
double GetFakeDataWeight_ND280ToNOvA_EnuKludge(int nu_pdg, int lep_pdg,
                                               int tgta, double E_nu_GeV,
                                               double PLep_GeV, double ThetaLep,
                                               int PrimSel,
                                               bool interpolate = true);

double GetFakeDataWeight_ND280ToNOvA_Enu(int nu_pdg, int lep_pdg, int tgta,
                                         double E_nu_GeV, int PrimSel,
                                         bool interpolate = true);

double GetFakeDataWeight_ND280ToNOvA_Q2(int nu_pdg, int lep_pdg, int tgta,
                                        double Q2_GeV2, int PrimSel,
                                        bool interpolate = true);
} // namespace t2knova

namespace t2knova {

enum reweightconfig {
  kT2KND_to_NOvA = 0,
  kT2KND_to_NOvA_EnuKludge,
  kT2KND_to_NOvA_Enu,
  kT2KND_to_NOvA_Q2,
  kNOvA_to_T2KND_plep,
  kNOvA_to_T2KND_Q2,
  kNOvA_to_T2KND_ptlep,
  kNoWeight
};

const char *all_rwconfig[] = {"t2knd_to_nova/EnuPLepThetaLep",
                              "DUMMY4KLUDGE",
                              "t2knd_to_nova/Enu",
                              "t2knd_to_nova/Q2",
                              "nova_to_t2k/EnuPLepEAvHad",
                              "nova_to_t2k/EnuQ2EAvHad",
                              "nova_to_t2k/EnuPtLepEAvHad"};

const char *all_tgta_str[] = {"H", "C", "O"};
const int all_tgta[] = {1, 12, 16};

static bool loaded = false;
static std::unordered_map<
    nuspecies,
    std::unordered_map<reweightconfig,
                       std::unordered_map<int, std::unique_ptr<TH1>>>>
    rwhists;
static std::unordered_map<
    nuspecies,
    std::unordered_map<reweightconfig,
                       std::unordered_map<int, std::unique_ptr<TH1>>>>
    EnuCorrections;

template <typename T, size_t N> inline size_t arrsize(T (&arr)[N]) { return N; }
template <typename T> inline size_t arrsize(T arr) { return arr.size(); }

inline void LoadHists(std::string const &inputfile = "FakeDataInputs.root") {
  std::unique_ptr<TFile> fin(new TFile(inputfile.c_str()));
  if (fin->IsZombie()) {
    std::cout << "Failed to read \"" << inputfile << "\"" << std::endl;
    abort();
  }
  int found = 0;
  for (int i = 0; i < arrsize(all_nuspecies); ++i) {
    nuspecies nuspec = nuspecies(i);
    std::string nuspec_str = all_nuspecies[i];

    for (int j = 0; j < arrsize(all_rwconfig); ++j) {
      reweightconfig rwconfig = reweightconfig(j);
      if (rwconfig ==
          kT2KND_to_NOvA_EnuKludge) { // this doesn't have its own histograms
        continue;
      }
      std::string rwconfig_str = all_rwconfig[j];

      for (selection sel : ReWeightSelectionList) {
        std::string sel_str = SelectionList[sel];

        for (int l = 0; l < arrsize(all_tgta_str); ++l) {
          std::string tgta_str = all_tgta_str[l];
          int tgta_sel_offset = all_tgta[l] * 100;

          rwhists[nuspec][rwconfig][tgta_sel_offset + sel] =
              GetTH1(fin, rwconfig_str + "/" + tgta_str + "/" + nuspec_str +
                              "/" + sel_str);
          if (rwhists[nuspec][rwconfig][tgta_sel_offset + sel]) {
            rwhists[nuspec][rwconfig][tgta_sel_offset + sel]->SetDirectory(
                nullptr);
            found++;
          } else {
            std::cout << "[WARN]: Failed to read "
                      << rwconfig_str + "/" + tgta_str + "/" + nuspec_str +
                             "/" + sel_str
                      << std::endl;
          }

          if (rwconfig == kT2KND_to_NOvA) {
            EnuCorrections[nuspec][rwconfig][tgta_sel_offset + sel] =
                GetTH1(fin, rwconfig_str + "/" + tgta_str + "/" + nuspec_str +
                                "/MissingPSENuCorrection_" + sel_str);
            if (EnuCorrections[nuspec][rwconfig][tgta_sel_offset + sel]) {
              EnuCorrections[nuspec][rwconfig][tgta_sel_offset + sel]
                  ->SetDirectory(nullptr);
            }
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

  std::unique_ptr<TH1> &rathist =
      rwhists[nuspec][kNOvA_to_T2KND_plep]
             [(tgta * 100) + (iscc ? kCCInc : kNCInc)];

  if (rathist) {
    double wght =
        EvalHist3D(rathist, E_nu_GeV, PLep_GeV, EVisHadronic_GeV, interpolate);
    return wght && !std::isnormal(wght) ? 0 : wght;
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

  std::unique_ptr<TH1> &rathist =
      rwhists[nuspec][kNOvA_to_T2KND_Q2]
             [(tgta * 100) + (iscc ? kCCInc : kNCInc)];
  if (rathist) {
    double wght =
        EvalHist3D(rathist, E_nu_GeV, Q2_GeV2, EVisHadronic_GeV, interpolate);
    return wght && !std::isnormal(wght) ? 0 : wght;
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

  std::unique_ptr<TH1> &rathist =
      rwhists[nuspec][kNOvA_to_T2KND_ptlep]
             [(tgta * 100) + (iscc ? kCCInc : kNCInc)];
  if (rathist) {
    double wght =
        EvalHist3D(rathist, E_nu_GeV, PtLep_GeV, EVisHadronic_GeV, interpolate);
    return wght && !std::isnormal(wght) ? 0 : wght;
  }
  return 1;
}

inline double GetFakeDataWeight_ND280ToNOvA(int nu_pdg, int lep_pdg, int tgta,
                                            double E_nu_GeV, double PLep_GeV,
                                            double ThetaLep, int PrimSel,
                                            bool interpolate) {
  if (!loaded) {
    LoadHists();
  }
  nuspecies nuspec = getnuspec(nu_pdg);

  bool iscc = (nu_pdg != lep_pdg);

  std::unique_ptr<TH1> &rathist =
      rwhists[nuspec][kT2KND_to_NOvA][(tgta * 100) + PrimSel];

  if (rathist) {
    return EvalHist3D(rathist, E_nu_GeV, PLep_GeV, ThetaLep, interpolate);
  }

  return 1;
}

inline double GetFakeDataWeight_ND280ToNOvA_EnuKludge(
    int nu_pdg, int lep_pdg, int tgta, double E_nu_GeV, double PLep_GeV,
    double ThetaLep, int PrimSel, bool interpolate) {
  if (!loaded) {
    LoadHists();
  }
  nuspecies nuspec = getnuspec(nu_pdg);

  bool iscc = (nu_pdg != lep_pdg);

  std::unique_ptr<TH1> &rathist =
      rwhists[nuspec][kT2KND_to_NOvA][(tgta * 100) + PrimSel];

  if (rathist) {
    return EvalHist3D(rathist, E_nu_GeV, PLep_GeV, ThetaLep, interpolate) *
           EvalHist1D(
               EnuCorrections[nuspec][kT2KND_to_NOvA][(tgta * 100) + PrimSel],
               E_nu_GeV, interpolate);
  }

  return 1;
}

inline double GetFakeDataWeight_ND280ToNOvA_Enu(int nu_pdg, int lep_pdg,
                                                int tgta, double E_nu_GeV,
                                                int PrimSel, bool interpolate) {
  if (!loaded) {
    LoadHists();
  }
  nuspecies nuspec = getnuspec(nu_pdg);

  bool iscc = (nu_pdg != lep_pdg);

  std::unique_ptr<TH1> &rathist =
      rwhists[nuspec][kT2KND_to_NOvA_Enu][(tgta * 100) + PrimSel];
  if (rathist) {
    return EvalHist1D(rathist, E_nu_GeV, interpolate);
  }

  return 1;
}

inline double GetFakeDataWeight_ND280ToNOvA_Q2(int nu_pdg, int lep_pdg,
                                               int tgta, double Q2_GeV2,
                                               int PrimSel, bool interpolate) {
  if (!loaded) {
    LoadHists();
  }
  nuspecies nuspec = getnuspec(nu_pdg);

  bool iscc = (nu_pdg != lep_pdg);

  std::unique_ptr<TH1> &rathist =
      rwhists[nuspec][kT2KND_to_NOvA_Q2][(tgta * 100) + PrimSel];
  if (rathist) {
    return EvalHist1D(rathist, Q2_GeV2, interpolate);
  }

  return 1;
}

} // namespace t2knova
