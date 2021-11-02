#include "ChannelHistCollections.h"
#include "T2KNOvATrueSelectionHelper.hxx"
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

using namespace t2knova;

bool bymode = false;

TH1F *totxsecs;

SelectionHists<TH1F> *Enu;
SelectionHists<TH1F> *PLep;
SelectionHists<TH1F> *ThetaLep;
SelectionHists<TH1F> *EAvHad;
SelectionHists<TH1F> *PtLep;
SelectionHists<TH1F> *Q2;
SelectionHists<TH1F> *q0;
SelectionHists<TH1F> *q3;

void Fill(TTreeReader &ttrdr, toml::value const &plots_config,
          t2knova::reweightconfig weightconfig, int tgta_select = 0) {
  Enu = SelectionHistsFromTOML<TH1F>("Enu", toml::find(plots_config, "Enu"));

  PLep = SelectionHistsFromTOML<TH1F>("PLep", toml::find(plots_config, "PLep"));
  Q2 = SelectionHistsFromTOML<TH1F>("Q2", toml::find(plots_config, "Q2"));
  q0 = SelectionHistsFromTOML<TH1F>("q0", toml::find(plots_config, "q0"));
  q3 = SelectionHistsFromTOML<TH1F>("q3", toml::find(plots_config, "q3"));

  ThetaLep = SelectionHistsFromTOML<TH1F>("ThetaLep",
                                          toml::find(plots_config, "ThetaLep"));

  EAvHad = SelectionHistsFromTOML<TH1F>("EAvHad",
                                        toml::find(plots_config, "EAvHad"));

  PtLep =
      SelectionHistsFromTOML<TH1F>("PtLep", toml::find(plots_config, "PtLep"));

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

  while (ttrdr.Next()) {
    if (ent_it && !(ent_it % shout_it)) {
      std::cout << "[Read] " << ent_it << "/" << nents << "("
                << (100 * ent_it / nents) << "%)" << std::endl;
    }

    if (tgta_select && (rdr.tgta() != tgta_select)) {
      ent_it++;
      continue;
    }

    double w = rdr.fScaleFactor() * rdr.RWWeight();

    int primary_selection = rdr.GetPrimarySelection();

    if (weightconfig == t2knova::kT2KND_to_NOvA) {
      w *= t2knova::GetFakeDataWeight_ND280ToNOvA(
          rdr.PDGNu(), rdr.PDGLep(), rdr.tgta(), rdr.Enu_true(), rdr.PLep(),
          acos(rdr.CosLep()), primary_selection);
    } else if (weightconfig == t2knova::kT2KND_to_NOvA_Enu) {
      w *= t2knova::GetFakeDataWeight_ND280ToNOvA_Enu(
          rdr.PDGNu(), rdr.PDGLep(), rdr.tgta(), rdr.Enu_true(),
          primary_selection);
    } else if (weightconfig == t2knova::kT2KND_to_NOvA_Q2) {
      w *= t2knova::GetFakeDataWeight_ND280ToNOvA_Q2(
          rdr.PDGNu(), rdr.PDGLep(), rdr.tgta(), rdr.Q2(), primary_selection);
    } else if (weightconfig == t2knova::kNOvA_to_T2KND_plep) {
      w *= t2knova::GetFakeDataWeight_NOvAToT2K_PLep(rdr.PDGNu(), rdr.PDGLep(),
                                                     rdr.tgta(), rdr.Enu_true(),
                                                     rdr.PLep(), rdr.EavAlt());
    } else if (weightconfig == t2knova::kNOvA_to_T2KND_Q2) {
      w *= t2knova::GetFakeDataWeight_NOvAToT2K_Q2(rdr.PDGNu(), rdr.PDGLep(),
                                                   rdr.tgta(), rdr.Enu_true(),
                                                   rdr.Q2(), rdr.EavAlt());
    } else if (weightconfig == t2knova::kNOvA_to_T2KND_ptlep) {
      w *= t2knova::GetFakeDataWeight_NOvAToT2K_PtLep(
          rdr.PDGNu(), rdr.PDGLep(), rdr.tgta(), rdr.Enu_true(),
          (rdr.PLep()) * sqrt(1 - pow(rdr.CosLep(), 2)), rdr.EavAlt());
    }

    std::vector<int> sels = rdr.GetSelections();

    if (!sels.size()) {
      ent_it++;
      continue;
    }

    Enu->Fill(w, sels, rdr.Mode(), rdr.Enu_true());
    PLep->Fill(w, sels, rdr.Mode(), rdr.PLep());
    ThetaLep->Fill(w, sels, rdr.Mode(), acos(rdr.CosLep()));
    EAvHad->Fill(w, sels, rdr.Mode(), rdr.EavAlt());
    PtLep->Fill(w, sels, rdr.Mode(),
                rdr.PLep() * sqrt(1 - pow(rdr.CosLep(), 2)));
    Q2->Fill(w, sels, rdr.Mode(), rdr.Q2());
    q0->Fill(w, sels, rdr.Mode(), rdr.q0());
    q3->Fill(w, sels, rdr.Mode(), rdr.q3());

    ent_it++;
  }
}

std::string input_file, output_file, output_dir;
std::string hist_config_file;

t2knova::reweightconfig wconfig = t2knova::kNoWeight;
int tgta_select = 0;

void SayUsage(char const *argv[]) {
  std::cout << "[USAGE]: " << argv[0]
            << "\n"
               "\t-i <inp.root>          : Input file"
               "\t-H <config.toml>       : RWHist config file"
               "\t-o <out.root>          : Output file"
               "\t-M                     : Separate by Mode"
               "\t-a <[C|H|O|any]>       : Target descriptor"
               "\t-W <Config>            : ReWeight Config"
               "\t         Configs:"
               "\t              * T2KND_to_NOvA_Enu"
               "\t              * T2KND_to_NOvA_Q2"
               "\t              * NOvA_to_T2KND_plep"
               "\t              * NOvA_to_T2KND_Q2"
               "\t              * NOvA_to_T2KND_ptlep"
               "\t-d </sub/dir/to/use>   : Output sub directory"
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
    } else if (std::string(argv[opt]) == "-H") {
      hist_config_file = argv[++opt];
    } else if (std::string(argv[opt]) == "-o") {
      output_file = argv[++opt];
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
                  << ". Should be ND280/NOvAND." << std::endl;
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
  gStyle->SetOptStat(false);

  handleOpts(argc, argv);

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
