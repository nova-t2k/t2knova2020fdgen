#define T2KNOVARW_MERGED_CC0PI

#include "ChannelHistCollections.h"

#include "T2KNOvATruthTreeReader.h"
#include "gstReader.h"

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "THn.h"
#include "TStyle.h"

#include <cmath>
#include <variant>

using namespace t2knova;

bool bymode = false;
bool isgst = false;

SelectionHists<TH1D> *Enu;
SelectionHists<TH1D> *PLep;
SelectionHists<TH1D> *q0;
SelectionHists<TH1D> *q3;

void Fill(TTreeReader &ttrdr, toml::value const &plots_config,
          int tgta_select = 0) {

  Enu = SelectionHistsFromTOML<TH1D>("Enu", plots_config);
  PLep = SelectionHistsFromTOML<TH1D>("PLep", plots_config);
  q0 = SelectionHistsFromTOML<TH1D>("q0", plots_config);
  q3 = SelectionHistsFromTOML<TH1D>("q3", plots_config);

  std::unique_ptr<T2KNOvATruthTreeReader> rdr_nuis(
      isgst ? nullptr : new T2KNOvATruthTreeReader(ttrdr));
  std::unique_ptr<gstReader> rdr_gst(isgst ? new gstReader(ttrdr) : nullptr);

  size_t nents = ttrdr.GetEntries(true);
  size_t ent_it = 0;
  size_t shout_it = nents / 25;

  double EnuCut = toml::find<double>(plots_config, "EnuCut");

  while (ttrdr.Next()) {

    int mode = bymode ? (isgst ? rdr_gst->Mode() : rdr_nuis->Mode()) : 0;

    if (ent_it && !(ent_it % shout_it)) {
      std::cout << "[Read] " << ent_it << "/" << nents << "("
                << (100 * ent_it / nents) << "%)" << std::endl;
    }

    if (tgta_select &&
        ((isgst ? rdr_gst->tgta() : rdr_nuis->tgta()) != tgta_select)) {
      ent_it++;
      continue;
    }

    std::vector<int> sels =
        (isgst ? rdr_gst->GetSelections() : rdr_nuis->GetSelections());

    if (!sels.size()) {
      ent_it++;
      continue;
    }

    int primary_selection = (isgst ? rdr_gst->GetPrimarySelection()
                                   : rdr_nuis->GetPrimarySelection());

    double Enu_gev = (isgst ? rdr_gst->Enu_true() : rdr_nuis->Enu_true());
    Enu->Fill(1.0, sels, mode, Enu_gev);
    PLep->Fill(1.0, sels, mode, (isgst ? rdr_gst->PLep() : rdr_nuis->PLep()));

    TLorentzVector fslep4 = (isgst ? rdr_gst->FSLepP4() : rdr_nuis->FSLepP4());
    TLorentzVector isnup4 = TLorentzVector(0, 0, Enu_gev, Enu_gev);
    TLorentzVector fourmomtransfer = (isnup4 - fslep4);

    q0->Fill(1.0, sels, mode, fourmomtransfer.E());
    q3->Fill(1.0, sels, mode, fourmomtransfer.Vect().Mag());

    ent_it++;
  }
}

std::string input_file = "", output_file = "", output_dir = "";
std::string hist_config_file = "";

int tgta_select = 0;

void SayUsage(char const *argv[]) {
  std::cout << "[USAGE]: " << argv[0]
            << "\n"
               "\n\t-i <inp.root>            : Input file"
               "\n\t--gst                    : Parse as GST file."
               "\n\t-H <config.toml>         : RWHist config file"
               "\n\t-o <out.root>            : Output file"
               "\n\t-M                       : Separate by Mode"
               "\n\t-a <[C|H|O|any]>         : Target descriptor"
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
        std::cout << "Invalid target selector passed: " << arg
                  << ". Should be C/H/O/any." << std::endl;
        SayUsage(argv);
        exit(1);
      }
    } else if (std::string(argv[opt]) == "-d") {
      output_dir = argv[++opt];
    } else if (std::string(argv[opt]) == "--gst") {
      isgst = true;
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

  TTreeReader ttrdr(isgst ? "gst" : "T2KNOvATruthTree", &fin);

  Fill(ttrdr, plots_config, tgta_select);

  TFile fout(output_file.c_str(), "UPDATE");
  if (fin.IsZombie()) {
    std::cout << "Failed to write " << output_file.c_str() << std::endl;
    return 2;
  }

  TDirectory *dout = MakeDirectoryStructure(&fout, output_dir);

  Enu->Write(dout, true);
  PLep->Write(dout, true);
  q0->Write(dout, true);
  q3->Write(dout, true);

  fout.Close();
}
