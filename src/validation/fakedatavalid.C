#include "ChannelHistCollections.h"
#include "T2KNOvA/TrueSelectionHelper.hxx"
#include "T2KNOvATruthTreeReader.h"
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

#include "OscillationHelper.hxx"

#include <cmath>

using namespace t2knova;

bool bymode = false;
bool dotune = true;
bool doosc = false;

TH1F *totxsecs;

SelectionHists<TH3F> *EnuPLepThetaLep;
SelectionHists<TH3F> *EnuPtLepEAvHad;
SelectionHists<TH1F> *Enu;
SelectionHists<TH1F> *ERecQE;
SelectionHists<TH1F> *PLep;
SelectionHists<TH1F> *ThetaLep;
SelectionHists<TH1F> *EAvHad;
SelectionHists<TH1F> *PtLep;
SelectionHists<TH1F> *Q2;
SelectionHists<TH1F> *q0;
SelectionHists<TH1F> *q3;
SelectionHists<TH2F> *q0q3;
SelectionHists<TH1F> *hmfscpip;
SelectionHists<TH1F> *hmfspi0p;
SelectionHists<TH1F> *ncpi;
SelectionHists<TH1F> *npi0;

OscillationHelper oh_disp, oh_app, oh_dispb, oh_appb;

double const mass_proton = 0.938272;
double const mass_neutron = 0.939565;

double EnuQErec(double elep, double plep, double costh, double binding,
                bool neutrino) {

  const double V = binding;       // binding potential
  const double mn = mass_neutron; // neutron mass
  const double mp = mass_proton;  // proton mass

  double mN_eff = mn - V;
  double mN_oth = mp;

  if (!neutrino) {
    mN_eff = mp - V;
    mN_oth = mn;
  }

  double el = elep;
  double pl = plep;                         // momentum of lepton
  double ml = std::sqrt(el * el - pl * pl); // lepton mass

  return (2 * mN_eff * el - ml * ml + mN_oth * mN_oth - mN_eff * mN_eff) /
         (2 * (mN_eff - el + pl * costh));
};

void Fill(TTreeReader &ttrdr, toml::value const &plots_config,
          t2knova::reweightconfig weightconfig, int tgta_select = 0) {
  EnuPLepThetaLep = SelectionHistsFromTOML<TH3F>(
      "EnuPLepThetaLep", toml::find(plots_config, "EnuPLepThetaLep"));
  EnuPtLepEAvHad = SelectionHistsFromTOML<TH3F>(
      "EnuPtLepEAvHad", toml::find(plots_config, "EnuPtLepEAvHad"));
  Enu = SelectionHistsFromTOML<TH1F>("Enu", toml::find(plots_config, "Enu"));

  ERecQE = SelectionHistsFromTOML<TH1F>("ERecQE", toml::find(plots_config, "ERecQE"));
  PLep = SelectionHistsFromTOML<TH1F>("PLep", toml::find(plots_config, "PLep"));
  Q2 = SelectionHistsFromTOML<TH1F>("Q2", toml::find(plots_config, "Q2"));
  q0 = SelectionHistsFromTOML<TH1F>("q0", toml::find(plots_config, "q0"));
  q3 = SelectionHistsFromTOML<TH1F>("q3", toml::find(plots_config, "q3"));
  q0q3 = SelectionHistsFromTOML<TH2F>("q0q3", toml::find(plots_config, "q0q3"));

  ThetaLep = SelectionHistsFromTOML<TH1F>("ThetaLep",
                                          toml::find(plots_config, "ThetaLep"));

  EAvHad = SelectionHistsFromTOML<TH1F>("EAvHad",
                                        toml::find(plots_config, "EAvHad"));

  PtLep =
      SelectionHistsFromTOML<TH1F>("PtLep", toml::find(plots_config, "PtLep"));

  hmfscpip = SelectionHistsFromTOML<TH1F>("hmfscpip",
                                          toml::find(plots_config, "hmfscpip"));
  hmfspi0p = SelectionHistsFromTOML<TH1F>("hmfspi0p",
                                          toml::find(plots_config, "hmfspi0p"));
  ncpi = SelectionHistsFromTOML<TH1F>("ncpi", toml::find(plots_config, "ncpi"));
  npi0 = SelectionHistsFromTOML<TH1F>("npi0", toml::find(plots_config, "npi0"));

  totxsecs = new TH1F("totxsecs", ";;#sigma^{#int#Phi} cm^{2}",
                      SelectionList.size(), 0, SelectionList.size());
  totxsecs->SetDirectory(nullptr);

  T2KNOvATruthTreeReader rdr(ttrdr);
  if (!bymode) {
    rdr.SetNoModeInfo();
  }

  size_t nents = ttrdr.GetEntries(true);
  size_t ent_it = 0;
  size_t shout_it = nents / 25;

  double EnuCut = toml::find<double>(plots_config, "EnuCut");

  while (ttrdr.Next()) {
    if (ent_it && !(ent_it % shout_it)) {
      std::cout << "[Read] " << ent_it << "/" << nents << "("
                << (100 * ent_it / nents) << "%)" << std::endl;
    }

    if (tgta_select && (rdr.tgta() != tgta_select)) {
      ent_it++;
      continue;
    }

    if (rdr.Enu_true() > EnuCut) {
      continue;
    }

    double w = rdr.fScaleFactor() * (dotune ? rdr.RWWeight() : 1);

    if (doosc && (std::abs(rdr.Mode()) < 30)) { //Oscillate CC events if enabled
      switch (rdr.PDGNu()) {
      case 14: {
        w *= oh_disp.GetWeight(rdr.Enu_true());
        break;
      }
      case -14: {
        w *= oh_dispb.GetWeight(rdr.Enu_true());
        break;
      }
      case 12: {
        w *= oh_app.GetWeight(rdr.Enu_true());
        break;
      }
      case -12: {
        w *= oh_appb.GetWeight(rdr.Enu_true());
        break;
      }
      default: {
      }
      }
    }

    int primary_selection = rdr.GetPrimarySelection();

    if (weightconfig == t2knova::kT2KND_to_NOvA) {
      w *= t2knova::GetFakeDataWeight_ND280ToNOvA(
          rdr.PDGNu(), rdr.PDGLep(), rdr.tgta(), rdr.Enu_true(), rdr.PLep(),
          rdr.AngLep_deg(), primary_selection, false);
    } else if (weightconfig == t2knova::kT2KND_to_NOvA_EnuKludge) {
      w *= t2knova::GetFakeDataWeight_ND280ToNOvA_EnuKludge(
          rdr.PDGNu(), rdr.PDGLep(), rdr.tgta(), rdr.Enu_true(), rdr.PLep(),
          rdr.AngLep_deg(), primary_selection, false);
    } else if (weightconfig == t2knova::kT2KND_to_NOvA_Enu) {
      w *= t2knova::GetFakeDataWeight_ND280ToNOvA_Enu(
          rdr.PDGNu(), rdr.PDGLep(), rdr.tgta(), rdr.Enu_true(),
          primary_selection, false);
    } else if (weightconfig == t2knova::kT2KND_to_NOvA_Q2) {
      w *= t2knova::GetFakeDataWeight_ND280ToNOvA_Q2(rdr.PDGNu(), rdr.PDGLep(),
                                                     rdr.tgta(), rdr.Q2(),
                                                     primary_selection, false);
    } else if (weightconfig == t2knova::kNOvA_to_T2KND_plep) {
      w *= t2knova::GetFakeDataWeight_NOvAToT2K_PLep(
          rdr.PDGNu(), rdr.PDGLep(), rdr.tgta(), rdr.Enu_true(), rdr.PLep(),
          rdr.Eav_NOvA(), primary_selection, false);
    } else if (weightconfig == t2knova::kNOvA_to_T2KND_Q2) {
      w *= t2knova::GetFakeDataWeight_NOvAToT2K_Q2(
          rdr.PDGNu(), rdr.PDGLep(), rdr.tgta(), rdr.Enu_true(), rdr.Q2(),
          rdr.Eav_NOvA(), primary_selection, false);
    } else if (weightconfig == t2knova::kNOvA_to_T2KND_ptlep) {
      w *= t2knova::GetFakeDataWeight_NOvAToT2K_PtLep(
          rdr.PDGNu(), rdr.PDGLep(), rdr.tgta(), rdr.Enu_true(),
          (rdr.PLep()) * sqrt(1 - pow(rdr.CosLep(), 2)), rdr.Eav_NOvA(),
          primary_selection, false);
    }

    std::vector<int> sels = rdr.GetSelections();

    if (!sels.size()) {
      ent_it++;
      continue;
    }

    EnuPLepThetaLep->Fill(w, sels, rdr.Mode(), rdr.Enu_true(), rdr.PLep(),
                          rdr.AngLep_deg());
    EnuPtLepEAvHad->Fill(w, sels, rdr.Mode(), rdr.Enu_true(),
                         rdr.PLep() * sqrt(1 - pow(rdr.CosLep(), 2)),
                         rdr.Eav_NOvA());
    Enu->Fill(w, sels, rdr.Mode(), rdr.Enu_true());
    ERecQE->Fill(w, sels, rdr.Mode(),
                 EnuQErec(rdr.FSLepP4().E(), rdr.PLep(), rdr.CosLep(), 0,
                          rdr.PDGNu() > 0));

    PLep->Fill(w, sels, rdr.Mode(), rdr.PLep());
    ThetaLep->Fill(w, sels, rdr.Mode(), rdr.AngLep_deg());
    EAvHad->Fill(w, sels, rdr.Mode(), rdr.Eav_NOvA());
    PtLep->Fill(w, sels, rdr.Mode(),
                rdr.PLep() * sqrt(1 - pow(rdr.CosLep(), 2)));
    Q2->Fill(w, sels, rdr.Mode(), rdr.Q2());
    q0->Fill(w, sels, rdr.Mode(), rdr.q0());
    q3->Fill(w, sels, rdr.Mode(), rdr.q3());
    q0q3->Fill(w, sels, rdr.Mode(), rdr.q0(), rdr.q3());
    hmfscpip->Fill(w, sels, rdr.Mode(), rdr.hmfscpip());
    hmfspi0p->Fill(w, sels, rdr.Mode(), rdr.hmfspi0p());
    ncpi->Fill(w, sels, rdr.Mode(), rdr.ncpi());
    npi0->Fill(w, sels, rdr.Mode(), rdr.npi0());

    ent_it++;
  }
}

std::string input_file, output_file, output_dir;
std::string hist_config_file;
std::string FDSInputs;

t2knova::reweightconfig wconfig = t2knova::kNoWeight;
int tgta_select = 0;

void SayUsage(char const *argv[]) {
  std::cout << "[USAGE]: " << argv[0]
            << "\n"
               "\n\t-i <inp.root>          : Input file"
               "\n\t-F <FDSInputs.root>    : FDS Weighting file"
               "\n\t-H <config.toml>       : RWHist config file"
               "\n\t-o <out.root>          : Output file"
               "\n\t-M                     : Separate by Mode"
               "\n\t--oscillate            : Apply oscillation weights"
               "\n\t--No-Tune              : Don't apply tune weights"
               "\n\t-a <[C|H|O|any]>       : Target descriptor"
               "\n\t-W <Config>            : ReWeight Config"
               "\n\t         Configs:"
               "\n\t              * T2KND_to_NOvA"
               "\n\t              * T2KND_to_NOvA_Enu"
               "\n\t              * T2KND_to_NOvA_Q2"
               "\n\t              * NOvA_to_T2KND_plep"
               "\n\t              * NOvA_to_T2KND_Q2"
               "\n\t              * NOvA_to_T2KND_ptlep"
               "\n\t-d </sub/dir/to/use>   : Output sub directory"
               " input tree. \n"
            << std::endl;
}

void handleOpts(int argc, char const *argv[]) {
  int opt = 1;
  while (opt < argc) {
    if ((std::string(argv[opt]) == "-?") ||
        (std::string(argv[opt]) == "--help")) {
      SayUsage(argv);
      exit(0);
    } else if (std::string(argv[opt]) == "-i") {
      input_file = argv[++opt];
    } else if (std::string(argv[opt]) == "-F") {
      FDSInputs = argv[++opt];
    } else if (std::string(argv[opt]) == "-H") {
      hist_config_file = argv[++opt];
    } else if (std::string(argv[opt]) == "-o") {
      output_file = argv[++opt];
    } else if (std::string(argv[opt]) == "--No-Tune") {
      dotune = false;
    } else if (std::string(argv[opt]) == "--oscillate") {
      doosc = true;
    } else if (std::string(argv[opt]) == "-a") {
      std::string arg = std::string(argv[++opt]);
      if (arg == "C") {
        tgta_select = 12;
      } else if (arg == "H") {
        tgta_select = 1;
      } else if (arg == "O") {
        tgta_select = 16;
      } else if (arg == "any") {
        tgta_select = 0;
      } else {
        std::cout << "Invalid target selector passed: " << arg
                  << ". Should be C/H/O/any." << std::endl;
        SayUsage(argv);
        exit(1);
      }
    } else if (std::string(argv[opt]) == "-W") {
      std::string arg = std::string(argv[++opt]);
      if (arg == "T2KND_to_NOvA") {
        wconfig = t2knova::kT2KND_to_NOvA;
      } else if (arg == "T2KND_to_NOvA_EnuKludge") {
        wconfig = t2knova::kT2KND_to_NOvA_EnuKludge;
      } else if (arg == "T2KND_to_NOvA_Enu") {
        wconfig = t2knova::kT2KND_to_NOvA_Enu;
      } else if (arg == "T2KND_to_NOvA_Q2") {
        wconfig = t2knova::kT2KND_to_NOvA_Q2;
      } else if (arg == "NOvA_to_T2KND_plep") {
        wconfig = t2knova::kNOvA_to_T2KND_plep;
      } else if (arg == "NOvA_to_T2KND_Q2") {
        wconfig = t2knova::kNOvA_to_T2KND_Q2;
      } else if (arg == "NOvA_to_T2KND_ptlep") {
        wconfig = t2knova::kNOvA_to_T2KND_ptlep;
      } else {
        std::cout << "Invalid weight config selector passed: " << arg
                  << std::endl;
        SayUsage(argv);
        exit(1);
      }
    } else if (std::string(argv[opt]) == "-d") {
      output_dir = argv[++opt];
    } else if (std::string(argv[opt]) == "-M") {
      bymode = true;
    }
    opt++;
  }
}

int main(int argc, char const *argv[]) {

  // Sin^2(Theta_12)
  double s2th12 = 0.297;
  // Sin^2(Theta_13)
  double s2th13 = 0.0214;
  // Sin^2(Theta_23)
  double s2th23 = 0.526;
  // Dm^2_21
  double dm2_21 = 7.37E-5;
  //|Dm^2_Atm|
  double dm2_atm = 2.463E-3;
  // dcp
  double dcp = 0;
  double osc_params[] = {s2th12, s2th13, s2th23, dm2_21, dm2_atm, dcp};
  oh_disp.Setup_baseline(osc_params, 295);
  oh_app.Setup_baseline(osc_params, 295);
  oh_disp.SetOscillationChannel(14, 14);
  oh_app.SetOscillationChannel(14, 12);
  oh_dispb.Setup_baseline(osc_params, 295);
  oh_appb.Setup_baseline(osc_params, 295);
  oh_dispb.SetOscillationChannel(-14, -14);
  oh_appb.SetOscillationChannel(-14, -12);

  gStyle->SetOptStat(false);

  handleOpts(argc, argv);

  t2knova::LoadHists(FDSInputs);

  toml::value const &plots_config_file = toml_h::parse_card(hist_config_file);
  toml::value const &plots_config =
      toml::find(plots_config_file, "FakeDataValidConfig");

  TFile fin(input_file.c_str());
  if (fin.IsZombie()) {
    std::cout << "Failed to read " << input_file << std::endl;
    return 2;
  }

  TTreeReader ttrdr("T2KNOvATruthTree", &fin);

  Fill(ttrdr, plots_config, wconfig, tgta_select);

  TFile fout(output_file.c_str(), "UPDATE");
  if (fin.IsZombie()) {
    std::cout << "Failed to write " << output_file.c_str() << std::endl;
    return 2;
  }

  TDirectory *dout = MakeDirectoryStructure(&fout, output_dir);

  for (int i = 0; i < SelectionList.size(); ++i) {
    totxsecs->GetXaxis()->SetBinLabel(i + 1, SelectionList[i].c_str());
  }

  dout->WriteTObject(totxsecs, "totxsecs");
  EnuPLepThetaLep->Write(dout, true);
  EnuPtLepEAvHad->Write(dout, true);
  Enu->Write(dout, true);
  PLep->Write(dout, true);
  ThetaLep->Write(dout, true);
  EAvHad->Write(dout, true);
  PtLep->Write(dout, true);
  Q2->Write(dout, true);
  q0->Write(dout, true);
  q3->Write(dout, true);
  q0q3->Write(dout, true);

  hmfscpip->Write(dout, true);
  hmfspi0p->Write(dout, true);
  ncpi->Write(dout, true);
  npi0->Write(dout, true);

  fout.Close();
}
