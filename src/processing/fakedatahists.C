#include "ChannelHistCollections.h"

#include "T2KNOvA/ROOTHelper.hxx"
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

bool bymode = false;

using namespace t2knova;

SelectionHists<TH3F> *EnuPLepThetaLep;
SelectionHists<TH3F> *EnuPtLepEAvHad;
TrueChannelHist<TH1F> *XSecs;

void Fill(TTreeReader &ttrdr, toml::value const &plots_config, bool ist2k,
          int tgta_select = 0) {

  if (ist2k) {
    EnuPLepThetaLep =
        SelectionHistsFromTOML<TH3F>("EnuPLepThetaLep", plots_config);
  } else {
    EnuPtLepEAvHad =
        SelectionHistsFromTOML<TH3F>("EnuPtLepEAvHad", plots_config);
  }

  XSecs = new TrueChannelHist<TH1F>("SelectionXSecs", ";Selection;Rate",
                                    AllSelectionList.size(), 0,
                                    AllSelectionList.size());

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

    std::vector<int> sels = rdr.GetSelections();

    if (!sels.size()) {
      ent_it++;
      continue;
    }

    if (ist2k) {
      EnuPLepThetaLep->Fill(w, sels, rdr.Mode(), rdr.Enu_true(), rdr.PLep(),
                            rdr.AngLep_deg());
    } else {
      EnuPtLepEAvHad->Fill(w, sels, rdr.Mode(), rdr.Enu_true(),
                           rdr.PLep() * sqrt(1 - pow(rdr.CosLep(), 2)),
                           rdr.Eav_NOvA());
    }

    for (auto sel : sels) {
      XSecs->Fill(w, rdr.Mode(), sel);
    }

    ent_it++;
  }
}

std::string input_file, output_file, output_dir;
std::string hist_config_file;
bool ist2k = false;
int tgta_select = 0;

void SayUsage(char const *argv[]) {
  std::cout << "[USAGE]: " << argv[0]
            << "\n"
               "\t-i <inp.root>          : Input file"
               "\t-H <config.toml>       : RWHist config file"
               "\t-o <out.root>          : Output file"
               "\t-M                     : Separate by Mode"
               "\t-a <[C|H|O|any]>       : Target descriptor"
               "\t-e <[ND280|NOvAND]>    : Experiment"
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
        std::cout << "Invalid target selector passed: " << argv[4]
                  << ". Should be C/H/O/any." << std::endl;
        abort();
      }
    } else if (std::string(argv[opt]) == "-e") {
      std::string arg = std::string(argv[++opt]);
      if (arg == "ND280") {
        ist2k = true;
      } else if (arg == "NOvAND") {
        ist2k = false;
      } else {
        std::cout << "Invalid experiment selector passed: " << argv[4]
                  << ". Should be ND280/NOvAND." << std::endl;
        abort();
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
      toml::find(plots_config_file, "FakeDataConfig");

  TFile fin(input_file.c_str());
  if (fin.IsZombie()) {
    std::cout << "Failed to read " << input_file << std::endl;
    return 2;
  }

  TTreeReader ttrdr("T2KNOvATruthTree", &fin);

  Fill(ttrdr, plots_config, ist2k, tgta_select);

  TFile fout(output_file.c_str(), "UPDATE");
  if (fin.IsZombie()) {
    std::cout << "Failed to write " << output_file.c_str() << std::endl;
    return 2;
  }

  TDirectory *dout = MakeDirectoryStructure(&fout, output_dir);

  XSecs->Apply([](TH1F &h) {
    for (int i = 0; i < SelectionList.size(); ++i) {
      h.GetXaxis()->SetBinLabel(i + 1, SelectionList[i].c_str());
    }
  });

  if (ist2k) {
    EnuPLepThetaLep->Write(dout, true);
  } else {
    EnuPtLepEAvHad->Write(dout, true);
  }

  XSecs->Write(dout);

  fout.Close();
}