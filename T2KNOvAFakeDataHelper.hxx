#pragma once

#include "TDirectory.h"
#include "TFile.h"
#include "TH3.h"

#include <map>

namespace t2knova {

TH1 *GetTH1(TFile *f, std::string const &name) {
  TDirectory *odir = gDirectory;

  TH1 *h;
  f->GetObject(name.c_str(), h);
  if (!h) {
    std::cout << "[ERROR]: Failed to find histogram named: " << name
              << std::endl;
  } else {
    h->SetDirectory(nullptr);
  }

  if (odir) {
    gDirectory->cd();
  }

  return h;
}

enum nuspecies { kNuMu, kNuMub, kNuE, kNuEb };
const char *all_nuspecies[] = {"numu", "numub", "nue", "nueb"};

nuspecies getnuspec(int pdg) {
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

enum reweightconfig { kT2KND_to_NOvA, kSK_to_NOvA, kNOvA_to_T2KND };
const char *all_rwconfig[] = {"t2knd_to_nova_", "sk_to_nova_", "nova_to_t2k_"};

enum selection {
  kCCINC,
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

bool loaded = false;
std::map<nuspecies, std::map<reweightconfig, std::map<selection, TH3D *>>>
    rwhists;

void LoadHists() {
  TFile fin(ifile.c_str());
  if (fin.IsZombie()) {
    std::cout << "Failed to read " << ifile.c_str() << std::endl;
    return 2;
  }

  for (int i = 0; i < 4; ++i) {
    std::string nuspec = all_nuspecies[i];

    for (int i = 0; i < 3; ++i) {
      std::string rwconfig = all_rwconfig[i];

      for (int i = 0; i < 10; ++i) {
        std::string sel = all_sel[i];

        rwhists[nuspec][rwconfig][sel] = dynamic_cast<TH3D *>(
            GetTH1(&fin, rwconfig + "_" + nuspec + "_" + sel));
        if (rwhists[nuspec][rwconfig][sel]) {
          rwhists[nuspec][rwconfig][sel]->SetDirectory(nullptr);
        }
      }
    }
  }
  loaded = true;
}

double GetFakeDataWeight_NOvAToT2K(int nu_pdg, int lep_pdg, double E_nu_GeV,
                                   double ELep_GeV, double EVisHadronic_GeV) {
  if (!loaded) {
    LoadHists();
  }
  nuspecies nuspec = getnuspec(nu_pdg);

  bool iscc = (nu_pdg == lep_pdg);

  TH3D *rathist = rwhists[nuspec][kNOvA_to_T2KND][iscc ? kCCINC : kNCINC];
  if (rathist) {

    int xbin = rathist->GetXaxis()->FindFixBin(E_nu_GeV);
    int ybin = rathist->GetYaxis()->FindFixBin(ELep_GeV);
    int zbin = rathist->GetZaxis()->FindFixBin(EVisHadronic_GeV);

    if ((xbin != 0) && (xbin != (rathist->GetXaxis()->GetNbins() + 1)) &&
        (ybin != 0) && (ybin != (rathist->GetYaxis()->GetNbins() + 1)) &&
        (zbin != 0) && (zbin != (rathist->GetZaxis()->GetNbins() + 1))) {
      return rathist->GetBinContent(xbin, ybin, zbin);
    }
  }

  return 1;
}

double GetFakeDataWeight_ND280ToNOvA(int nu_pdg, int lep_pdg, double E_nu_GeV,
                                     double ELep_GeV, double ThetaLep,
                                     int NFSCPi, int NFSPi0, int NOther) {
  if (!loaded) {
    LoadHists();
  }
  nuspecies nuspec = getnuspec(nu_pdg);

  bool iscc = (nu_pdg == lep_pdg);

  selection sel = kCCINC;
  if ((NFSCPi + NFSPi0 + NOther) == 0) {
    sel = iscc ? kCC0pi : kNC0pi;
  } else if ((NFSCPi + NFSPi0 + NOther) == 1) {
    if (NFSCPi == 1) {
      sel = iscc ? kCC1cpi : kNC1cpi;
    } else if (NFSCPi == 1) {
      sel = iscc ? kCC1pi0 : kNC1pi0;
    } else {
      sel = iscc ? kCCOther : kNCOther;
    }
  } else {
    sel = iscc ? kCCOther : kNCOther;
  }

  TH3D *rathist = rwhists[nuspec][kT2KND_to_NOvA][sel];
  if (rathist) {

    int xbin = rathist->GetXaxis()->FindFixBin(E_nu_GeV);
    int ybin = rathist->GetYaxis()->FindFixBin(ELep_GeV);
    int zbin = rathist->GetZaxis()->FindFixBin(ThetaLep);

    if ((xbin != 0) && (xbin != (rathist->GetXaxis()->GetNbins() + 1)) &&
        (ybin != 0) && (ybin != (rathist->GetYaxis()->GetNbins() + 1)) &&
        (zbin != 0) && (zbin != (rathist->GetZaxis()->GetNbins() + 1))) {
      return rathist->GetBinContent(xbin, ybin, zbin);
    }
  }

  return 1;
}

double GetFakeDataWeight_SKToNOvA(int nu_pdg, int lep_pdg, double E_nu_GeV,
                                  double ELep_GeV, double ThetaLep, int NFSCPi,
                                  int NFSPi0, int NOther) {
  if (!loaded) {
    LoadHists();
  }
  nuspecies nuspec = getnuspec(nu_pdg);

  bool iscc = (nu_pdg == lep_pdg);

  selection sel = kCCINC;
  if ((NFSCPi + NFSPi0 + NOther) == 0) {
    sel = iscc ? kCC0pi : kNC0pi;
  } else if ((NFSCPi + NFSPi0 + NOther) == 1) {
    if (NFSCPi == 1) {
      sel = iscc ? kCC1cpi : kNC1cpi;
    } else if (NFSCPi == 1) {
      sel = iscc ? kCC1pi0 : kNC1pi0;
    } else {
      sel = iscc ? kCCOther : kNCOther;
    }
  } else {
    sel = iscc ? kCCOther : kNCOther;
  }

  TH3D *rathist = rwhists[nuspec][kSK_to_NOvA][sel];
  if (rathist) {

    int xbin = rathist->GetXaxis()->FindFixBin(E_nu_GeV);
    int ybin = rathist->GetYaxis()->FindFixBin(ELep_GeV);
    int zbin = rathist->GetZaxis()->FindFixBin(ThetaLep);

    if ((xbin != 0) && (xbin != (rathist->GetXaxis()->GetNbins() + 1)) &&
        (ybin != 0) && (ybin != (rathist->GetYaxis()->GetNbins() + 1)) &&
        (zbin != 0) && (zbin != (rathist->GetZaxis()->GetNbins() + 1))) {
      return rathist->GetBinContent(xbin, ybin, zbin);
    }
  }

  return 1;
}

} // namespace t2knova