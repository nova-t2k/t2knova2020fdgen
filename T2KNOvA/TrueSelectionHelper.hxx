#pragma once

#include <algorithm>
#include <iostream>
#include <set>
#include <sstream>
#include <string>
#include <vector>

// #define DEBUG_FSP_LOOP

namespace t2knova {

struct FSParticleSummary {
  FSParticleSummary()
      : IsCC(false), PrimLepPDG(0), NChLep(0), NNeutralLep(0), NProton(0),
        NNeutron(0), NAntiNucleon(0), NChPi(0), NPi0(0), NGamma(0), EGamma(0),
        NChKaons(0), NNeutralKaons(0), NOtherMesons(0), NOtherBaryons(0),
        NOther(0), NNuclear(0) {}

  bool IsCC;
  int PrimLepPDG;
  int NChLep;
  int NNeutralLep;
  int NProton;
  int NNeutron;
  int NAntiNucleon;
  int NChPi;
  int NPi0;
  int NGamma;
  double EGamma;
  int NChKaons;
  int NNeutralKaons;
  int NOtherMesons;
  int NOtherBaryons;
  int NOther;
  int NNuclear;
  int NOvAFSIMode;

  int GetNExoticHadrons() {
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
    count += NChKaons;
    count += NNeutralKaons;
    count += NOtherMesons;
    count += NOtherBaryons;
    count += NOther;

    return count;
  }

  int GetNPions() { return NChPi + NPi0; }

  std::string Print() const {
    std::stringstream ss;

    ss << "{ IsCC: " << IsCC << ", PrimLepPDG: " << PrimLepPDG
       << ", NChLep: " << NChLep << ", NNeutralLep: " << NNeutralLep
       << ", NProton: " << NProton << ", NNeutron: " << NNeutron
       << ", NAntiNucleon: " << NAntiNucleon << ", NChPi: " << NChPi
       << ", NPi0: " << NPi0 << ", NGamma: " << NGamma
       << ", NChKaons: " << NChKaons << ", NNeutralKaons: " << NNeutralKaons
       << ", NOtherMesons: " << NOtherMesons
       << ", NOtherBaryons: " << NOtherBaryons << ", NOther: " << NOther
       << ", NNuclear: " << NNuclear << "}";
    return ss.str();
  }
};

static std::set<int> UnknownParticles;

template <typename T>
inline FSParticleSummary
T2KNOvAFlatTreeToFSParticleSummary(int NFSP, int *FSPDG, T *FSE_GeV,
                                   int NOvAFSIMode) {
  FSParticleSummary fsps;
  fsps.NOvAFSIMode = NOvAFSIMode;

#ifdef DEBUG_FSP_LOOP
  std::cout << "NFSP: " << NFSP << std::endl;
#endif

  for (int i = 0; i < NFSP; ++i) {
#ifdef DEBUG_FSP_LOOP
    std::cout << "\t" << i << ": PDG = " << FSPDG[i] << ", E = " << FSE_GeV[i]
              << std::endl;
#endif

    switch (FSPDG[i]) {
    case 11:
    case 13:
    case 15:
    case -11:
    case -13:
    case -15: {
      if (!fsps.PrimLepPDG) {
        fsps.PrimLepPDG = FSPDG[i];
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
      if (!fsps.PrimLepPDG) {
        fsps.PrimLepPDG = FSPDG[i];
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
      fsps.EGamma += FSE_GeV[i];
      break;
    }
    case 321:
    case 311:
    case -321:
    case -311: {
      fsps.NChKaons++;
      break;
    }
    case 310:
    case 130: {
      fsps.NNeutralKaons++;
      break;
    }

    case 221: {
      fsps.NOtherMesons++;
      break;
    }

    case 3212: /*Sigma0*/
    case 3222: /*Sigma+*/
    case 3112: /*Sigma-*/
    case -3122:
    case 3122: /*Lambda*/ {
      fsps.NOtherBaryons++;
      break;
    }

    default: {
      if (FSPDG[i] > 1000000000) {
        fsps.NNuclear++;
      } else {
        if (!UnknownParticles.count(FSPDG[i])) {
          std::cout << "[INFO]: Found Other particle: " << FSPDG[i]
                    << std::endl;
          UnknownParticles.insert(FSPDG[i]);
        }
        fsps.NOther++;
      }
      break;
    }
    }
  }
  fsps.IsCC = fsps.PrimLepPDG % 2;

  return fsps;
}

#define SEL_LIST                                                               \
  SEL_X(CCInc)                                                                 \
  SEL_X(CCInc_RW)                                                              \
  SEL_X(CC0pi)                                                                 \
  SEL_X(CC0pi_QE)                                                              \
  SEL_X(CC0pi_2p2h)                                                            \
  SEL_X(CC0pi_Other)                                                           \
  SEL_X(CC1Gamma)                                                              \
  SEL_X(CCDeExciteGamma)                                                       \
  SEL_X(CCNGamma)                                                              \
  SEL_X(CC1cpi)                                                                \
  SEL_X(CC1pi0)                                                                \
  SEL_X(CCmultipi)                                                             \
  SEL_X(CCOther)                                                               \
  SEL_X(CCOther_QE)                                                            \
  SEL_X(NCInc)                                                                 \
  SEL_X(NCInc_RW)                                                              \
  SEL_X(NC0pi)                                                                 \
  SEL_X(NC1Gamma)                                                              \
  SEL_X(NCDeExciteGamma)                                                       \
  SEL_X(NCNGamma)                                                              \
  SEL_X(NC1cpi)                                                                \
  SEL_X(NC1pi0)                                                                \
  SEL_X(NCmultipi)                                                             \
  SEL_X(NCOther)                                                               \
  SEL_X(NCOther_QE)                                                            \
  SEL_X(NuEElastic)                                                            \
  SEL_X(NOvAFSIMode_CC0Pi)                                                     \
  SEL_X(NOvAFSIMode_CC1Pi)                                                     \
  SEL_X(NOvAFSIMode_CCOth)                                                     \
  SEL_X(NOvAFSIMode_CCInc)                                                     \
  SEL_X(NOvAFSIMode_NCInc)                                                     \
  SEL_X(NOvAFSIMode_CC1cPi)                                                    \
  SEL_X(NOvAFSIMode_CC1Pi0)                                                    \
  SEL_X(NOvAFSIMode_CCMultiPi)                                                 \
  SEL_X(NOvAFSIMode_Nope)                                                      \
  SEL_X(NoPrimarySel)

#define X_CATNAME(A, B) A##B
#define SEL_X(a) X_CATNAME(k, a),

enum selection { SEL_LIST };
static std::vector<selection> AllSelectionList = {SEL_LIST};

#undef SEL_X
#define SEL_X(a) #a,

static std::vector<std::string> SelectionList = {SEL_LIST};

#undef X_CATNAME
#undef SEL_X

#undef SEL_LIST

#ifdef T2KNOVARW_MERGED_CC0PI
static std::vector<selection> ReWeightSelectionList = {
    kCC0pi, kCC1cpi, kCC1pi0, kCCmultipi, kCCOther,
    kNC0pi, kNC1cpi, kNC1pi0, kNCmultipi, kNCOther};
#else
static std::vector<selection> ReWeightSelectionList = {
    kCC0pi_QE, kCC0pi_2p2h, kCC0pi_Other, kCC1cpi, kCC1pi0,    kCCmultipi,
    kCCOther,  kNC0pi,      kNC1cpi,      kNC1pi0, kNCmultipi, kNCOther};
#endif

inline std::vector<int> GetSelections(FSParticleSummary fsps, int mode = 0) {
  std::vector<int> sels = std::vector<int>{fsps.IsCC ? kCCInc : kNCInc};

  sels.push_back(fsps.NOvAFSIMode);

  const double kEGammaDeExciteCut = 10 * 1E-3;

  // if EGamma < 10 MeV then we will call it a de-excitation
  if ((fsps.EGamma < kEGammaDeExciteCut) && (fsps.NGamma > 0)) {
    sels.push_back(fsps.IsCC ? kCCDeExciteGamma : kNCDeExciteGamma);
  }

  // Normalish event with just nucleons and pions
  if (fsps.GetNExoticHadrons() == 0) {

    // Filter out gamma-ie events
    if (fsps.EGamma > kEGammaDeExciteCut) {
      sels.push_back(fsps.IsCC ? ((fsps.NGamma > 1) ? kCCNGamma : kCC1Gamma)
                               : ((fsps.NGamma > 1) ? kNCNGamma : kNC1Gamma));
    } else {
      // If it doesn't have gammas then we're reweighting it
      sels.push_back(fsps.IsCC ? kCCInc_RW : kNCInc_RW);

      if (fsps.GetNPions() == 0) {

        // Have to do this first so that it is considered the primary selection
        // if the mode info is available
        if (mode && fsps.IsCC) {
          int CC0pi_offset = std::abs(mode);
          if (std::abs(mode) > 2) {
            CC0pi_offset = 3;
          }
          sels.push_back(kCC0pi + CC0pi_offset);
        }
        sels.push_back(fsps.IsCC ? kCC0pi : kNC0pi);

      } else if ((fsps.NPi0 == 0) && ((fsps.NChPi) == 1)) {
        sels.push_back(fsps.IsCC ? kCC1cpi : kNC1cpi);
      } else if ((fsps.NPi0 == 1) && ((fsps.NChPi) == 0)) {
        sels.push_back(fsps.IsCC ? kCC1pi0 : kNC1pi0);
      } else {
        sels.push_back(fsps.IsCC ? kCCmultipi : kNCmultipi);
      }
    }
    return sels;
  } else if (fsps.GetNExoticHadrons() == 1) {
    int NNeutralLep_extra = fsps.NNeutralLep - (fsps.IsCC ? 0 : 1);
    if ((NNeutralLep_extra == 1) && (std::abs(fsps.PrimLepPDG) == 11)) {
      return {
          kNuEElastic,
      };
    }
  }

  if (std::abs(mode) == 1) {
    sels.push_back(kCCOther_QE);
  } else {
    sels.push_back(fsps.IsCC ? kCCInc_RW : kNCInc_RW);
    sels.push_back(fsps.IsCC ? kCCOther : kNCOther);
  }

  return sels;
}

inline int GetPrimarySelection(FSParticleSummary fsps, int Mode = 0) {
  auto sels = GetSelections(fsps, Mode);
  for (auto s : sels) {
    if (std::find(ReWeightSelectionList.begin(), ReWeightSelectionList.end(),
                  s) != ReWeightSelectionList.end()) {
      return s;
    }
  }
  return kNoPrimarySel;
}

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

} // namespace t2knova