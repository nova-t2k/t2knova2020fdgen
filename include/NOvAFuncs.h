#pragma once

#include "T2KNOvA/TrueSelectionHelper.hxx"

#include "TLorentzVector.h"

namespace t2knova {

template <typename T>
struct vect4 {
  vect4(std::array<T, 4> const &vec) {
    TLorentzVector v(vec[0], vec[1], vec[2], vec[3]);
    E = v.E();
    _gamma = v.Gamma();
  }
  T E;
  T _gamma;
  T Gamma() { return _gamma; }
};

template <typename T>
struct nova_part {
  nova_part(int _pdg, std::array<T, 4> const &vec) : pdg(_pdg), p(vec) {}
  int pdg;
  vect4<T> p;
};

template <typename T>
selection NOvAFSIMode(int nfsp, int *pdg, T *px, T *py, T *pz, T *E) {
  bool first_lep = true;

  int nfsipip = 0;
  int nfsipim = 0;
  int nfsipi0 = 0;
  int nfsilep = 0;
  int nfsiphn = 0;
  int nfsinuc = 0;
  int nfsioth = 0;
  for (int i = 0; i < nfsp; ++i) {
    nova_part<T> particle(pdg[i], std::array<T, 4>{px[i], py[i], pz[i], E[i]});

    if (particle.pdg > 1000000000) {
      return kNOvAFSIMode_Nope; // COH interaction off of nuclei
    }

    if (particle.pdg == 211) {
      nfsipip++;
    } else if (particle.pdg == -211) {
      nfsipim++;
    } else if (particle.pdg == 111) {
      nfsipi0++;
    } else if (particle.pdg == 22) {
      nfsiphn++;
    } else if ((abs(particle.pdg) > 10 && abs(particle.pdg) < 17)) {
      if (first_lep && !(abs(particle.pdg) & 1)) {
        return kNOvAFSIMode_NCInc;
      }
      first_lep = false;
      nfsilep++;
    } else if (particle.pdg == 2112 || particle.pdg == 2212) {
      nfsinuc++;
    } else {
      nfsioth++;
    }
  }

  if (nfsilep > 1) {
    nfsioth += nfsilep - 1;
  }

  if (nfsipip + nfsipim + nfsipi0 + nfsioth == 0) {
    return kNOvAFSIMode_CC0Pi;
  } else if (nfsipip + nfsipim == 1 && nfsipi0 + nfsioth == 0) {
    return kNOvAFSIMode_CC1cPi;
  } else if (nfsipi0 == 1 && nfsipip + nfsipim + nfsioth == 0) {
    return kNOvAFSIMode_CC1Pi0;
  } else if (nfsipip + nfsipim + nfsipi0 > 1 && nfsioth == 0) {
    return kNOvAFSIMode_CCMultiPi;
  } else {
    return kNOvAFSIMode_CCOth;
  }
}

template <typename T>
T Eav_NOvA(int nfsp, int *pdg, T *px, T *py, T *pz, T *E) {

  T eAvail = 0;

  for (int i = 0; i < nfsp; ++i) {
    nova_part<T> particle(pdg[i], std::array<T, 4>{px[i], py[i], pz[i], E[i]});

    // Code below here is from Zoya (unaltered from NOvA/CAFAna
    // implementation)

    double particle_Eavail = 0;
    // protons, pions: KE only
    if (particle.pdg == 2212 || abs(particle.pdg) == 211) {
      double gamma = particle.p.Gamma();
      particle_Eavail = (gamma - 1) / gamma * particle.p.E;
    }
    // pi0s, electrons, photons: total energy
    else if (particle.pdg == 111 || abs(particle.pdg) == 11 ||
             particle.pdg == 22)
      particle_Eavail = particle.p.E;

    else if (particle.pdg >= 2000000000) {
      // skip the bindinos
    }

    else if (particle.pdg >= 1000000000) {
      // do nothing for nucleons
    } else if (particle.pdg >= 2000 && particle.pdg != 2212 &&
               particle.pdg != 2112) {
      particle_Eavail = particle.p.E - 0.9382;
      // Primarily strange baryons add total energy minus proton mass since
      // decays will mostly contain protons
    } else if (particle.pdg <= -2000) {
      particle_Eavail = particle.p.E + 0.9382;
      // Primarily anti-protons add total energy plus proton mass since
      // anhillation is mostly the interaction mode
    } else if (particle.pdg != 2112 &&
               (abs(particle.pdg) < 11 ||
                abs(particle.pdg) > 16)) { // no neutrons or leptons
      particle_Eavail = particle.p.E;
      // mostly kaons add all the energy
    }

    eAvail += particle_Eavail;
  } // end loop over primaries

  return eAvail;
}

} // namespace t2knova