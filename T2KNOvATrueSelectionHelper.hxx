#pragma once

#include <iostream>
#include <set>
#include <sstream>
#include <string>
#include <vector>

// #define DEBUG_FSP_LOOP

namespace t2knova {

struct FSParticleSummary {
  FSParticleSummary()
      : IsCC(false),
        NChLep(0),
        NNeutralLep(0),
        NProton(0),
        NNeutron(0),
        NAntiNucleon(0),
        NChPi(0),
        NPi0(0),
        NGamma(0),
        NChKaons(0),
        NNeutralKaons(0),
        NOther(0),
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
  int NChKaons;
  int NNeutralKaons;
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
    count += NChKaons;
    count += NNeutralKaons;
    count += NOther;

    return count;
  }

  std::string Print() const {
    std::stringstream ss;

    ss << "{ IsCC: " << IsCC << ", NChLep: " << NChLep
       << ", NNeutralLep: " << NNeutralLep << ", NProton: " << NProton
       << ", NNeutron: " << NNeutron << ", NAntiNucleon: " << NAntiNucleon
       << ", NChPi: " << NChPi << ", NPi0: " << NPi0 << ", NGamma: " << NGamma
       << ", NChKaons: " << NChKaons << ", NNeutralKaons: " << NNeutralKaons
       << ", NOther: " << NOther << ", NNuclear: " << NNuclear << "}";
    return ss.str();
  }
};

static std::set<int> UnknownParticles;

inline FSParticleSummary T2KNOvAFlatTreeToFSParticleSummary(int NFSP,
                                                            int *FSPDG) {
  int PDG_FSLep = 0;

  FSParticleSummary fsps;

#ifdef DEBUG_FSP_LOOP
  std::cout << "NFSP: " << NFSP << std::endl;
#endif

  for (int i = 0; i < NFSP; ++i) {
#ifdef DEBUG_FSP_LOOP
    std::cout << "\t" << i << " " << FSPDG[i] << std::endl;
#endif

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
        fsps.NChKaons++;
        break;
      }
      case 310:
      case 130: {
        fsps.NNeutralKaons++;
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
  fsps.IsCC = PDG_FSLep % 2;

  return fsps;
}

#define SEL_LIST   \
  SEL_X(CCInc)     \
  SEL_X(CC0pi)     \
  SEL_X(CC1Gamma)  \
  SEL_X(CCNGamma)  \
  SEL_X(CC1cpi)    \
  SEL_X(CC1pi0)    \
  SEL_X(CCmultipi) \
  SEL_X(CCOther)   \
  SEL_X(NCInc)     \
  SEL_X(NC0pi)     \
  SEL_X(NC1Gamma)  \
  SEL_X(NCNGamma)  \
  SEL_X(NC1cpi)    \
  SEL_X(NC1pi0)    \
  SEL_X(NCmultipi) \
  SEL_X(NCOther)

#define X_CATNAME(A, B) A##B
#define SEL_X(a) X_CATNAME(k, a),

enum selection { SEL_LIST };

#undef SEL_X
#define SEL_X(a) #a,

std::vector<std::string> SelectionList = {SEL_LIST};

#undef SEL_LIST
#undef X_CATNAME
#undef SEL_X

inline std::vector<int> GetSelections(FSParticleSummary fsps) {
  std::vector<int> sels;

  sels.push_back(fsps.IsCC ? kCCInc : kNCInc);

  if (fsps.GetNNonNucleonPions() ==
      0) {  // Normalish event with just nucleons and pions
    if ((fsps.NPi0 + fsps.NChPi) == 0) {
      sels.push_back(fsps.IsCC ? kCC0pi : kNC0pi);
    } else if ((fsps.NPi0 == 0) && ((fsps.NChPi) == 1)) {
      sels.push_back(fsps.IsCC ? kCC1cpi : kNC1cpi);
    } else if ((fsps.NPi0 == 1) && ((fsps.NChPi) == 0)) {
      sels.push_back(fsps.IsCC ? kCC1pi0 : kNC1pi0);
    } else {
      sels.push_back(fsps.IsCC ? kCCmultipi : kNCmultipi);
    }
  } else if (fsps.GetNNonNucleonPions() == fsps.NGamma) {
    if (fsps.NGamma == 1) {
      sels.push_back(fsps.IsCC ? kCC1Gamma : kNC1Gamma);
    } else {
      sels.push_back(fsps.IsCC ? kCCNGamma : kNCNGamma);
    }
  } else {
    sels.push_back(fsps.IsCC ? kCCOther : kNCOther);
  }

#ifdef DEBUG_FSP_LOOP
  std::cout << fsps.Print() << " is selected as: " << std::endl;
  for (int s : sels) {
    std::cout << "\t" << SelectionList[s] << "(" << s << ")" << std::endl;
  }
#endif

  return sels;
}

inline int GetPrimarySelection(FSParticleSummary fsps) {
  auto sels = GetSelections(fsps);
  for (auto s : sels) {
    if ((s != kCCInc) && (s != kNCInc)) {
      return s;
    }
  }
  throw;
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

}  // namespace t2knova