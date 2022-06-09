#pragma once

#include <iostream>
#include <sstream>
#include <unordered_map>
#include <cmath>

#include "T2KNOvA/ROOTHelper.hxx"
#include "T2KNOvA/TrueSelectionHelper.hxx"

#include "TDirectory.h"
#include "TFile.h"
#include "TH3.h"
#include "THn.h"

// #define DEBUG_HIST_LOADER

namespace t2knova {

constexpr double NMinEvs = 1;
const double MaxFracError = 1.0 / std::sqrt(NMinEvs);

double GetFakeDataWeight_NOvAToT2K_PtLep(int nu_pdg, int lep_pdg, int tgta,
                                         double E_nu_GeV, double PtLep_GeV,
                                         double EVisHadronic_GeV,
                                         int PrimSel = -1,
                                         bool interpolate = true);
double GetFakeDataWeight_ND280ToNOvA(int nu_pdg, int lep_pdg, int tgta,
                                     double E_nu_GeV, double PLep_GeV,
                                     double ThetaLep, int PrimSel,
                                     bool interpolate = true);
} // namespace t2knova

namespace t2knova {

enum reweightconfig { kT2KND_to_NOvA = 0, kNOvA_to_T2KND_ptlep, kNoWeight };

const char *all_rwconfig[] = {"t2knd_to_nova/EnuPLepThetaLep",
                              "nova_to_t2k/EnuPtLepEAvHad"};

const char *all_tgta_str[] = {"H", "C", "O"};
const int all_tgta[] = {1, 12, 16};

static bool loaded = false;
static std::unordered_map<
    int, std::unordered_map<int, std::unordered_map<int, std::unique_ptr<TH1>>>>
    rwhists;

template <typename T, size_t N> inline size_t arrsize(T (&arr)[N]) { return N; }
template <typename T> inline size_t arrsize(T arr) { return arr.size(); }

inline void LoadHists(std::string const &inputfile = "FakeDataInputs.root") {
  std::unique_ptr<TFile> fin(new TFile(inputfile.c_str()));
  if (fin->IsZombie()) {
    std::cout << "Failed to read \"" << inputfile << "\"" << std::endl;
    abort();
  }
  int found = 0;
  for (size_t i = 0; i < arrsize(all_nuspecies); ++i) {
    nuspecies nuspec = nuspecies(i);
    std::string nuspec_str = all_nuspecies[i];

    for (size_t j = 0; j < arrsize(all_rwconfig); ++j) {
      reweightconfig rwconfig = reweightconfig(j);

      std::string rwconfig_str = all_rwconfig[j];

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

inline double GetFakeDataWeight_NOvAToT2K_PtLep(int nu_pdg, int lep_pdg,
                                                int tgta, double E_nu_GeV,
                                                double PtLep_GeV,
                                                double EVisHadronic_GeV,
                                                int PrimSel, bool interpolate) {
  if (!loaded) {
    LoadHists();
  }

  if (PrimSel == kNoPrimarySel) {
    return 1;
  }

  nuspecies nuspec = getnuspec(nu_pdg);

  bool iscc = (nu_pdg != lep_pdg);

  int offset = (PrimSel == -1) ? (iscc ? kCCInc : kNCInc) : PrimSel;

  std::unique_ptr<TH1> &rathist =
      rwhists[nuspec][kNOvA_to_T2KND_ptlep][(tgta * 100) + offset];

  // std::cout << "[RW]: iscc: " << iscc << ", offset: " << offset
  //           << ", PrimSel: " << PrimSel << std::endl;
  // std::cout << "[RW] rathist: " << rathist.get() << ", E_nu_GeV: " << E_nu_GeV
  //           << ", PtLep_GeV: " << PtLep_GeV
  //           << ", EVisHadronic_GeV: " << EVisHadronic_GeV << std::endl;

  if (rathist) {
    double wght =
        EvalHist3D(rathist, E_nu_GeV, PtLep_GeV, EVisHadronic_GeV, interpolate);
    // std::cout << "--wght: " << wght << std::endl;
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

  if (PrimSel == kNoPrimarySel) {
    return 1;
  }

  nuspecies nuspec = getnuspec(nu_pdg);

  std::unique_ptr<TH1> &rathist =
      rwhists[nuspec][kT2KND_to_NOvA][(tgta * 100) + PrimSel];

  if (rathist) {
    return EvalHist3D(rathist, E_nu_GeV, PLep_GeV, ThetaLep, interpolate);
  }

  return 1;
}

} // namespace t2knova
