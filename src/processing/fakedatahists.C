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

bool HasSel(std::vector<int> const &sels, int sel_to_check) {
  return (std::find(sels.begin(), sels.end(), sel_to_check) != sels.end());
}

enum FDS { kGenerated, kNDTuned, kMnv1Pi, kNonQE };

SelectionHists<TH3D> *EnuPLepThetaLep;
SelectionHists<TH3D> *EnuPtLepEAvHad;
TrueChannelHist<TH1D> *XSecs;
FDS FDSSet = kGenerated;

void Fill(TTreeReader &ttrdr, toml::value const &plots_config, bool ist2k,
          int tgta_select = 0) {

  if (ist2k) {
    EnuPLepThetaLep =
        SelectionHistsFromTOML<TH3D>("EnuPLepThetaLep", plots_config);
  } else {
    EnuPtLepEAvHad =
        SelectionHistsFromTOML<TH3D>("EnuPtLepEAvHad", plots_config);
  }

  XSecs = new TrueChannelHist<TH1D>("SelectionXSecs", ";Selection;Rate",
                                    AllSelectionList.size(), 0,
                                    AllSelectionList.size());

  T2KNOvATruthTreeReader rdr(ttrdr);

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

    double w = rdr.fScaleFactor();

    std::vector<int> sels = rdr.GetSelections();

    if (!sels.size()) {
      ent_it++;
      continue;
    }

    if (FDSSet != kGenerated) {
      w *= rdr.RWWeight();

      if (FDSSet == kMnv1Pi) {
        if ((std::abs(rdr.Mode()) >= 11) && (std::abs(rdr.Mode()) <= 13)) {
          w *= GetMINERvASPPLowQ2SuppressionWeight(rdr.Q2());
        }
      } else if (FDSSet == kNonQE) {
        if (std::abs(rdr.Mode()) == 1) {
          w *= UnWeightQ2BinWeights_T2K2020(rdr.Q2());
        } else if (std::find(sels.begin(), sels.end(), kCC0pi) != sels.end()) {
          w *= GetnonQEWeight(rdr.PDGNu(), rdr.GetQ2QE());
        }
      }
    }

    int mode = bymode ? rdr.Mode() : 0;

    if (ist2k) {
      EnuPLepThetaLep->Fill(w, sels, mode, rdr.Enu_true(), rdr.PLep(),
                            rdr.AngLep_deg());
    } else {
      EnuPtLepEAvHad->Fill(w, sels, mode, rdr.Enu_true(),
                           rdr.PLep() * sqrt(1 - pow(rdr.CosLep(), 2)),
                           rdr.Eav_NOvA());
    }

    for (auto sel : sels) {
      XSecs->Fill(w, bymode ? rdr.Mode() : 0, sel);
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
               "\n\t-i <inp.root>          : Input file"
               "\n\t-H <config.toml>       : RWHist config file"
               "\n\t-o <out.root>          : Output file"
               "\n\t--FDS <None|Mnv1Pi|NonQE> : Additional FDS weights to apply"
               "\n\t-M                     : Separate by Mode"
               "\n\t-a <[C|H|O|any]>       : Target descriptor"
               "\n\t-e <[ND280|NOvAND]>    : Experiment"
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
    } else if (std::string(argv[opt]) == "--FDS") {
      std::string arg = std::string(argv[++opt]);
      if (arg == "Generated") {
        FDSSet = kGenerated;
      } else if (arg == "NDTuned") {
        FDSSet = kNDTuned;
      } else if (arg == "Mnv1Pi") {
        FDSSet = kMnv1Pi;
      } else if (arg == "NonQE") {
        FDSSet = kNonQE;
      } else {
        std::cout << "Invalid FDS selector passed: " << argv[4]
                  << ". Should be Generated/NDTuned/Mnv1Pi/NonQE." << std::endl;
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

  XSecs->Apply([](TH1D &h) {
    for (int i = 0; i < SelectionList.size(); ++i) {
      h.GetXaxis()->SetBinLabel(i + 1, SelectionList[i].c_str());
    }
  });

  if (ist2k) {
    EnuPLepThetaLep->Write(dout, true);
    auto Enu = EnuPLepThetaLep->Transform<TH1D>(
        [](TH3D const &in) -> TH1D { return TH1D(*in.ProjectionX()); });
    Enu.SetName("Enu");
    Enu.Write(dout, true);
  } else {
    EnuPtLepEAvHad->Write(dout, true);
    auto Enu = EnuPtLepEAvHad->Transform<TH1D>(
        [](TH3D const &in) -> TH1D { return TH1D(*in.ProjectionX()); });
    Enu.SetName("Enu");
    Enu.Write(dout, true);
  }

  XSecs->Write(dout);

  fout.Close();
}