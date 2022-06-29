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

enum FDS { kGenerated, kNDTuned, kMnv1Pi, kNonQE };

///\note See https://arxiv.org/abs/1903.01558
inline double GetMINERvASPPLowQ2SuppressionWeight(double Q2_True_GeV,
                                                  double parameter_value = 1) {

  // Fit parameters for FrInel+Low-Q2 tune were
  // MA_RES = 0.93 +/- 0.05
  // NormRes = 116 +/- 7
  // NonRes1pi = 46 +/- 4
  // NonRes2pi = 120 +/- 32
  // ThetaPi = 1.0 (at limit)
  // FrInel = 132 +/- 27
  // R1 = 0.37 +/- 0.09
  // R2 = 0.60 +/- 0.16
  //
  // Fit parameters for FrAbs+Low-Q2 tune
  // MA_RES = 0.92 +/- 0.02
  // NormRes = 116 +/- 3
  // NonRes1pi = 46 +/- 4
  // NonRes2pi = 99 +/- 31
  // ThetaPi = 1.0 (at limit)
  // FrAbs = 48 +/- 21
  // R1 = 0.32 +/- 0.06
  // R2 = 0.5 (limit)
  static double const Q2_Max = 0.7;
  static double const Q2_t1 = 0;
  static double const Q2_t2 = 0.35;
  static double const R1 = 0.37;
  static double const R2 = 0.6;

  if ((Q2_True_GeV > Q2_Max) || (Q2_True_GeV < 0)) {
    return 1;
  }

  double RQ2 = (R2 * ((Q2_True_GeV - Q2_t1) * (Q2_True_GeV - Q2_Max)) /
                ((Q2_t2 - Q2_t1) * (Q2_t2 - Q2_Max))) +
               (((Q2_True_GeV - Q2_t1) * (Q2_True_GeV - Q2_t2)) /
                ((Q2_Max - Q2_t1) * (Q2_Max - Q2_t2)));
  return 1 - parameter_value * ((1 - R1) * pow((1 - RQ2), 2));
}

// BANFFEventBase::
inline double GetTrueEnuQE(T2KNOvATruthTreeReader &rdr) {
  static double const proton_mass = 0.938272;  // GeV
  static double const neutron_mass = 0.939566; // GeV
  bool nu = (rdr.PDGNu() > 0);
  double Eb = 27. * 0.001; // GeV
  double m_in = neutron_mass - Eb * 1.e-3;
  double m_out = proton_mass;
  if (!nu) {
    m_in = proton_mass - Eb * 1.e-3;
    m_out = neutron_mass;
  }
  auto FSLepP4 = rdr.FSLepP4();
  double muon_energy = FSLepP4.E();
  return (2. * m_in * muon_energy - FSLepP4.Mag2() - m_in * m_in +
          m_out * m_out) /
         (2. * (m_in - muon_energy +
                FSLepP4.Vect().Mag() * FSLepP4.Vect().CosTheta()));
}

// BANFFEventBase::
double GetQ2QE(T2KNOvATruthTreeReader &rdr) {
  auto FSLepP4 = rdr.FSLepP4();
  double muon_energy = FSLepP4.E();
  return -FSLepP4.Mag2() +
         2. * GetTrueEnuQE(rdr) *
             (muon_energy - FSLepP4.Vect().Mag() * FSLepP4.Vect().CosTheta());
}

inline double UnWeightQ2BinWeights_T2K2020(double Q2_True_GeV) {
  static double Q2Weights[] = {0.7841, 0.8868, 1.0228, 1.0268,
                               1.0867, 1.2568, 1.1360, 1.2593};
  static double fBinWidth_GeV2 = 0.05;
  double Q2Weight = 1;
  if (Q2_True_GeV <= 0.25) {
    Q2Weight = Q2Weights[int(std::floor(Q2_True_GeV / fBinWidth_GeV2))];
  } else if (Q2_True_GeV > 0.25 && Q2_True_GeV <= 0.5) {
    Q2Weight = Q2Weights[5];
  } else if (Q2_True_GeV > 0.5 && Q2_True_GeV <= 1.0) {
    Q2Weight = Q2Weights[6];
  } else if (Q2_True_GeV > 1.0) {
    Q2Weight = Q2Weights[7];
  }
  return 1.0 / Q2Weight;
}

inline double GetnonQEWeight(int nuPDG, double Q2_Reco_GeV) {

  static bool first = true;
  static TH1 *nuWeights = nullptr;
  static TH1 *nubWeights = nullptr;

  if (first) {
    char *T2KNOVA_INPUTS = getenv("T2KNOVA_INPUTS");

    if (!T2KNOVA_INPUTS) {
      std::cout << "[ERROR]: Expected T2KNOVA_INPUTS environment variable to "
                   "be defined."
                << std::endl;
      abort();
    }

    TFile *fin = new TFile(
        (std::string(T2KNOVA_INPUTS) + "/ScalingHisto_nu_antinu.root").c_str());
    if (!fin || !fin->IsOpen()) {
      std::cout << "[ERROR]; Failed to open " << T2KNOVA_INPUTS
                << "/ScalingHisto_nu_antinu.root" << std::endl;
      abort();
    }

    fin->GetObject("ScalingHisto_FGD1_numu", nuWeights);
    fin->GetObject("ScalingHisto_FGD1_anumu", nubWeights);

    nuWeights = dynamic_cast<TH1 *>(nuWeights->Clone("ScalingHisto_FGD1_numu"));
    if (!nubWeights) {
      std::cout << "[ERROR]: Failed to read ScalingHisto_FGD1_numu from "
                << T2KNOVA_INPUTS << "/ScalingHisto_nu_antinu.root"
                << std::endl;
      abort();
    }
    nuWeights->SetDirectory(nullptr);
    nubWeights =
        dynamic_cast<TH1 *>(nubWeights->Clone("ScalingHisto_FGD1_anumu"));
    if (!nubWeights) {
      std::cout << "[ERROR]: Failed to read ScalingHisto_FGD1_anumu from "
                << T2KNOVA_INPUTS << "/ScalingHisto_nu_antinu.root"
                << std::endl;
      abort();
    }
    nubWeights->SetDirectory(nullptr);

    fin->Close();
    delete fin;
    first = false;
  }

  TH1 *h = nuPDG > 0 ? nuWeights : nubWeights;

  int bin = h->FindFixBin(Q2_Reco_GeV);

  return h->GetBinContent(bin);
}

SelectionHists<TH3F> *EnuPLepThetaLep;
SelectionHists<TH3F> *EnuPtLepEAvHad;
TrueChannelHist<TH1F> *XSecs;
FDS FDSSet = kGenerated;

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
        if ((std::find(sels.begin(), sels.end(), kCC1cpi) != sels.end()) ||
            (std::find(sels.begin(), sels.end(), kCC1pi0) != sels.end())) {
          w *= GetMINERvASPPLowQ2SuppressionWeight(rdr.Q2());
        }
      } else if (FDSSet == kNonQE) {
        if (rdr.Mode() == 1) {
          w *= UnWeightQ2BinWeights_T2K2020(rdr.Q2());
        } else if (std::find(sels.begin(), sels.end(), kCC0pi) != sels.end()) {
          w *= GetnonQEWeight(rdr.PDGNu(), GetQ2QE(rdr));
        }
      }
    }

    if (ist2k) {
      EnuPLepThetaLep->Fill(w, sels, bymode ? rdr.Mode() : 0, rdr.Enu_true(),
                            rdr.PLep(), rdr.AngLep_deg());
    } else {
      EnuPtLepEAvHad->Fill(w, sels, bymode ? rdr.Mode() : 0, rdr.Enu_true(),
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
               "\t-i <inp.root>          : Input file"
               "\t-H <config.toml>       : RWHist config file"
               "\t-o <out.root>          : Output file"
               "\t--FDS <None|Mnv1Pi|NonQE> : Additional FDS weights to apply"
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
                  << ". Should be None/Mnv1Pi/NonQE." << std::endl;
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