#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "THn.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TProfile.h"
#include "TStyle.h"
#include "TTreeReader.h"

#include "T2KNOvAFakeDataHelper.hxx"

using namespace t2knova;

TH1D *totxsecs;

hblob<TH1D> *Enu;
hblob<TH1D> *PLep;
hblob<TH1D> *ThetaLep;
hblob<TH1D> *EAvHad;
hblob<TH1D> *PtLep;
hblob<TH1D> *Q2;
hblob<TH1D> *q0;
hblob<TH1D> *q3;

void Fill(TTreeReader &rdr,
          t2knova::reweightconfig weightconfig = t2knova::kNoWeight,
          int tgta_select = 0) {
  Enu = new hblob<TH1D>(
      "Enu",
      ";#it{E_{#nu}} (GeV);#it{y}; d#sigma/d#it{E_{#nu}} cm^{2} GeV^{-1}", 100,
      0, 5);

  PLep = new hblob<TH1D>(
      "PLep",
      ";#it{p}_{lep} (GeV/#it{c});#it{y}; d#sigma/d#it{p}_{lep} "
      "cm^{2} GeV^{-1} #it{c}",
      100, 0, 5);
  Q2 = new hblob<TH1D>("Q2",
                       ";#Q^{2} (GeV^{2});#it{y}; d#sigma/d#Q^{2} "
                       "cm^{2} GeV^{-2}",
                       100, 0, 5);
  q0 = new hblob<TH1D>("q0",
                       ";q_{0} (GeV);#it{y}; d#sigma/dq_{0} "
                       "cm^{2} GeV^{-1} #it{c}",
                       100, 0, 3);
  q3 = new hblob<TH1D>("q3",
                       ";q_{3} (GeV/#it{c});#it{y}; d#sigma/dq_{3} "
                       "cm^{2} GeV^{-1} #it{c}",
                       100, 0, 3);

  ThetaLep = new hblob<TH1D>(
      "ThetaLep",
      ";#theta_{lep} (GeV);#it{y}; d#sigma/d#theta_{lep} cm^{2} rad^{-1}", 100,
      0, M_PI);

  EAvHad = new hblob<TH1D>("EAvHad",
                           ";#it{E}_{had}^{vis} (GeV);#it{y}; "
                           "d#sigma/d#it{E}_{had}^{vis} cm^{2} GeV^{-1}",
                           100, 0, 5);

  PtLep = new hblob<TH1D>("PtLep",
                          ";#{p}^{T}_{lep} (GeV);#it{y}; "
                          "d#sigma/d#{p}^{T}_{lep} cm^{2} GeV^{-1} #it{c}",
                          100, 0, 2);

  totxsecs = new TH1D("totxsecs", ";;#sigma^{#int#Phi} cm^{2}", 10, 0, 10);
  totxsecs->SetDirectory(nullptr);

  TTreeReaderValue<int> PDGnu(rdr, "PDGnu");
  TTreeReaderValue<int> PDGLep(rdr, "PDGLep");
  TTreeReaderValue<int> Mode(rdr, "Mode");

  TTreeReaderValue<float> EavAlt(rdr, "EavAlt");
  TTreeReaderValue<float> Enu_true(rdr, "Enu_true");
  TTreeReaderValue<float> PLep_v(rdr, "PLep");
  TTreeReaderValue<float> CosLep(rdr, "CosLep");
  TTreeReaderValue<float> Q2_v(rdr, "Q2");
  TTreeReaderValue<float> q0_v(rdr, "q0");
  TTreeReaderValue<float> q3_v(rdr, "q3");
  TTreeReaderValue<int> tgta(rdr, "tgta");

  TTreeReaderValue<int> NFSpip(rdr, "NFSpip");
  TTreeReaderValue<int> NFSpim(rdr, "NFSpim");
  TTreeReaderValue<int> NFSpi0(rdr, "NFSpi0");
  TTreeReaderValue<int> NFSgamma(rdr, "NFSgamma");
  TTreeReaderValue<int> NFSlep(rdr, "NFSlep");
  TTreeReaderValue<int> NFSOther(rdr, "NFSOther");

  TTreeReaderValue<double> fScaleFactor(rdr, "fScaleFactor");
  TTreeReaderValue<double> RWWeight(rdr, "RWWeight");

  size_t nents = rdr.GetEntries(true);
  size_t ent_it = 0;
  size_t shout_it = nents / 25;

  while (rdr.Next()) {

    if (ent_it && !(ent_it % shout_it)) {
      std::cout << "[Read] " << ent_it << "/" << nents << "("
                << (100 * ent_it / nents) << "%)" << std::endl;
    }

    if (tgta_select && (*tgta != tgta_select)) {
      ent_it++;
      continue;
    }

    double w = *fScaleFactor * *RWWeight;

    int NCpi = *NFSpip + *NFSpim;
    int NOther = (*NFSlep > 0 ? (*NFSlep - 1) : 0) + *NFSOther;

    if (weightconfig == t2knova::kT2KND_to_NOvA) {
      w *= t2knova::GetFakeDataWeight_ND280ToNOvA(
          *PDGnu, *PDGLep, *tgta, *Enu_true, *PLep_v, acos(*CosLep), NCpi,
          *NFSpi0, NOther);
    } else if (weightconfig == t2knova::kT2KND_to_NOvA_Enu) {
      w *= t2knova::GetFakeDataWeight_ND280ToNOvA_Enu(
          *PDGnu, *PDGLep, *tgta, *Enu_true, NCpi, *NFSpi0, NOther);
    } else if (weightconfig == t2knova::kT2KND_to_NOvA_Q2) {
      w *= t2knova::GetFakeDataWeight_ND280ToNOvA_Q2(
          *PDGnu, *PDGLep, *tgta, *Q2_v, NCpi, *NFSpi0, NOther);
    } else if (weightconfig == t2knova::kNOvA_to_T2KND_plep) {
      w *= t2knova::GetFakeDataWeight_NOvAToT2K_PLep(
          *PDGnu, *PDGLep, *tgta, *Enu_true, *PLep_v, *EavAlt);
    } else if (weightconfig == t2knova::kNOvA_to_T2KND_Q2) {
      w *= t2knova::GetFakeDataWeight_NOvAToT2K_Q2(*PDGnu, *PDGLep, *tgta,
                                                   *Enu_true, *Q2_v, *EavAlt);
    } else if (weightconfig == t2knova::kNOvA_to_T2KND_ptlep) {
      w *= t2knova::GetFakeDataWeight_NOvAToT2K_PtLep(
          *PDGnu, *PDGLep, *tgta, *Enu_true,
          (*PLep_v) * sqrt(1 - pow(*CosLep, 2)), *EavAlt);
    }

    bool iscc = (*PDGLep) % 2;

    int NFSCPi = *NFSpip + *NFSpim;
    int NFSOther_gammaextralep = *NFSOther + *NFSgamma;
    if((*NFSlep) > 1){
      NFSOther_gammaextralep += (*NFSlep) - 1;
    }

    FlagBlob fblob = GetFlagBlob(iscc, NFSCPi, *NFSpi0, NFSOther_gammaextralep, *Mode);

    if (fblob.flagCCINC) {
      totxsecs->Fill(0.0, w);
      if (fblob.flagCC0pi) {
        totxsecs->Fill(1.0, w);
      } else if (fblob.flagCC1cpi) {
        totxsecs->Fill(2.0, w);
      } else if (fblob.flagCC1pi0) {
        totxsecs->Fill(3.0, w);
      } else {
        totxsecs->Fill(4.0, w);
      }
    } else if (fblob.flagNCINC) {
      totxsecs->Fill(5.0, w);
      if (fblob.flagNC0pi) {
        totxsecs->Fill(6.0, w);
      } else if (fblob.flagNC1cpi) {
        totxsecs->Fill(7.0, w);
      } else if (fblob.flagNC1pi0) {
        totxsecs->Fill(8.0, w);
      } else {
        totxsecs->Fill(9.0, w);
      }
    }

    Enu->Fill(w, fblob, *Enu_true);

    PLep->Fill(w, fblob, *PLep_v);
    ThetaLep->Fill(w, fblob, acos(*CosLep));
    EAvHad->Fill(w, fblob, *EavAlt);
    PtLep->Fill(w, fblob, *PLep_v * sqrt(1 - pow(*CosLep, 2)));
    Q2->Fill(w, fblob, *Q2_v);
    q0->Fill(w, fblob, *q0_v);
    q3->Fill(w, fblob, *q3_v);
    ent_it++;
  }
}

int main(int argc, char const *argv[]) {
  if (argc < 3) {
    std::cout << "Expects 2 arguments." << std::endl;
    return 1;
  }

  gStyle->SetOptStat(false);

  TFile fin(argv[1]);
  if (fin.IsZombie()) {
    std::cout << "Failed to read " << argv[1] << std::endl;
    return 2;
  }

  TTreeReader rdr("T2KNOvATruthTree", &fin);

  t2knova::reweightconfig weightconfig = t2knova::kNoWeight;

  if (std::string(argv[2]) == "NOvA_to_T2KND_plep") {
    weightconfig = t2knova::kNOvA_to_T2KND_plep;
  } else if (std::string(argv[2]) == "NOvA_to_T2KND_Q2") {
    weightconfig = t2knova::kNOvA_to_T2KND_Q2;
  } else if (std::string(argv[2]) == "NOvA_to_T2KND_ptlep") {
    weightconfig = t2knova::kNOvA_to_T2KND_ptlep;
  } else if (std::string(argv[2]) == "T2KND_To_NOvA") {
    weightconfig = t2knova::kT2KND_to_NOvA;
  } else if (std::string(argv[2]) == "T2KND_To_NOvA_Enu") {
    weightconfig = t2knova::kT2KND_to_NOvA_Enu;
  } else if (std::string(argv[2]) == "T2KND_To_NOvA_Q2") {
    weightconfig = t2knova::kT2KND_to_NOvA_Q2;
  }

  int tgta_select = 0;
  if (argc > 5) {
    if (std::string(argv[5]) == "C") {
      tgta_select = 12;
    } else if (std::string(argv[5]) == "H") {
      tgta_select = 1;
    } else if (std::string(argv[5]) == "O") {
      tgta_select = 16;
    } else if (std::string(argv[5]) == "any") {
      tgta_select = 0;
    } else {
      std::cout << "Invalid target selector passed: " << argv[4]
                << ". Should be C/H/O/any." << std::endl;
      return 1;
    }
  }

  Fill(rdr, weightconfig, tgta_select);

  TFile fout(argv[3], "UPDATE");
  if (fin.IsZombie()) {
    std::cout << "Failed to write " << argv[3] << std::endl;
    return 2;
  }

  TDirectory *dout = &fout;
  if (argc > 4) {
    std::string fqdir = argv[4];
    while (fqdir.find("/") != std::string::npos) {
      std::string ndir = fqdir.substr(0, fqdir.find("/"));
      if (!dout->GetDirectory(ndir.c_str())) {
        dout = dout->mkdir(ndir.c_str());
      } else {
        dout = dout->GetDirectory(ndir.c_str());
      }
      fqdir = fqdir.substr(fqdir.find("/") + 1);
    }
    dout = dout->mkdir(fqdir.c_str());
  }

  totxsecs->GetXaxis()->SetBinLabel(1, "CCInc");
  totxsecs->GetXaxis()->SetBinLabel(2, "CC0#pi");
  totxsecs->GetXaxis()->SetBinLabel(3, "CC1#pi^{#pm}");
  totxsecs->GetXaxis()->SetBinLabel(4, "CC1#pi^{0}");
  totxsecs->GetXaxis()->SetBinLabel(5, "CCMulti#pi + CCOther");
  totxsecs->GetXaxis()->SetBinLabel(6, "NCInc");
  totxsecs->GetXaxis()->SetBinLabel(7, "NC0#pi");
  totxsecs->GetXaxis()->SetBinLabel(8, "NC1#pi^{#pm}");
  totxsecs->GetXaxis()->SetBinLabel(9, "NC1#pi^{0}");
  totxsecs->GetXaxis()->SetBinLabel(10, "NCMulti#pi + CCOther");
  dout->WriteTObject(totxsecs, "totxsecs");

  Enu->Write(dout);
  PLep->Write(dout);
  ThetaLep->Write(dout);
  EAvHad->Write(dout);
  PtLep->Write(dout);
  Q2->Write(dout);
  q0->Write(dout);
  q3->Write(dout);

  fout.Close();
}
