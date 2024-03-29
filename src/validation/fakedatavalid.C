#define T2KNOVARW_MERGED_CC0PI

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

#include <cmath>

using namespace t2knova;

bool bymode = false;

enum FDS { kGenerated, kNDTuned, kMnv1Pi, kNonQE };
FDS FDSSet = kGenerated;

TrueChannelHist<TH1D> *XSecs;

SelectionHists<TH1D> *Enu;
SelectionHists<TH1D> *ERecQE;
SelectionHists<TH1D> *PLep;
SelectionHists<TH1D> *ThetaLep;
SelectionHists<TH1D> *CosThetaLep;
SelectionHists<TH2D> *PThetaLep;
SelectionHists<TH2D> *PThetaLep_outlier_low;
SelectionHists<TH2D> *PThetaLep_outlier_high;
SelectionHists<TH1D> *EAvHad;
SelectionHists<TH1D> *PtLep;
SelectionHists<TH1D> *Q2;
SelectionHists<TH1D> *q0;
SelectionHists<TH1D> *q3;
SelectionHists<TH1D> *yrec;
SelectionHists<TH2D> *Enuyrec;
SelectionHists<TH2D> *EnuQ2;
SelectionHists<TH2D> *EnuQ2_outlier_low;
SelectionHists<TH2D> *EnuQ2_outlier_high;
SelectionHists<TH2D> *EnuERecQE;
SelectionHists<TH2D> *EnuERecQEBias;
SelectionHists<TH2D> *EnuERecAvBias;
SelectionHists<TH2D> *q0_QEEAvBias;
SelectionHists<TH2D> *q0q3_high_outlier_low;
SelectionHists<TH2D> *q0q3_high_outlier_high;
SelectionHists<TH2D> *q0q3_low;
SelectionHists<TH2D> *q0q3_high;
SelectionHists<TH2D> *EnuEAvHad;
SelectionHists<TH1D> *hmfscpip;
SelectionHists<TH1D> *hmfspi0p;
SelectionHists<TH1D> *ncpi;
SelectionHists<TH1D> *npi0;
SelectionHists<TH1D> *hmfsprotonp;
SelectionHists<TH1D> *hmfsneutronp;
SelectionHists<TH1D> *nproton;
SelectionHists<TH1D> *nneutron;
SelectionHists<TH1D> *EGamma;
SelectionHists<TH1D> *EGamma_DeExcite;

void Fill(TTreeReader &ttrdr, toml::value const &plots_config,
          t2knova::reweightconfig weightconfig, int tgta_select = 0) {

  Enu = SelectionHistsFromTOML<TH1D>("Enu", plots_config);

  ERecQE = SelectionHistsFromTOML<TH1D>("ERecQE", plots_config);
  PLep = SelectionHistsFromTOML<TH1D>("PLep", plots_config);
  Q2 = SelectionHistsFromTOML<TH1D>("Q2", plots_config);
  EnuQ2 = SelectionHistsFromTOML<TH2D>("EnuQ2", plots_config);
  EnuERecQE = SelectionHistsFromTOML<TH2D>("EnuERecQE", plots_config);
  EnuERecQEBias = SelectionHistsFromTOML<TH2D>("EnuERecQEBias", plots_config);
  EnuERecAvBias = SelectionHistsFromTOML<TH2D>("EnuERecAvBias", plots_config);
  q0_QEEAvBias = SelectionHistsFromTOML<TH2D>("q0_QEEAvBias", plots_config);

  q0q3_low = SelectionHistsFromTOML<TH2D>("q0q3_low", plots_config);
  q0q3_high = SelectionHistsFromTOML<TH2D>("q0q3_high", plots_config);
  EnuEAvHad = SelectionHistsFromTOML<TH2D>("EnuEAvHad", plots_config);

  q0 = SelectionHistsFromTOML<TH1D>("q0", plots_config);
  q3 = SelectionHistsFromTOML<TH1D>("q3", plots_config);
  yrec = SelectionHistsFromTOML<TH1D>("yrec", plots_config);
  Enuyrec = SelectionHistsFromTOML<TH2D>("Enuyrec", plots_config);

  ThetaLep = SelectionHistsFromTOML<TH1D>("ThetaLep", plots_config);
  CosThetaLep = SelectionHistsFromTOML<TH1D>("CosThetaLep", plots_config);

  PThetaLep = SelectionHistsFromTOML<TH2D>("PThetaLep", plots_config);

  EAvHad = SelectionHistsFromTOML<TH1D>("EAvHad", plots_config);

  PtLep = SelectionHistsFromTOML<TH1D>("PtLep", plots_config);

  hmfscpip = SelectionHistsFromTOML<TH1D>("hmfscpip", plots_config);
  hmfspi0p = SelectionHistsFromTOML<TH1D>("hmfspi0p", plots_config);
  ncpi = SelectionHistsFromTOML<TH1D>("ncpi", plots_config);
  npi0 = SelectionHistsFromTOML<TH1D>("npi0", plots_config);

  hmfsprotonp = SelectionHistsFromTOML<TH1D>("hmfsprotonp", plots_config);
  hmfsneutronp = SelectionHistsFromTOML<TH1D>("hmfsneutronp", plots_config);
  nproton = SelectionHistsFromTOML<TH1D>("nproton", plots_config);
  nneutron = SelectionHistsFromTOML<TH1D>("nneutron", plots_config);

  EGamma = SelectionHistsFromTOML<TH1D>("EGamma", plots_config);
  EGamma_DeExcite =
      SelectionHistsFromTOML<TH1D>("EGamma_DeExcite", plots_config);

  XSecs = new TrueChannelHist<TH1D>("SelectionXSecs", ";Selection;Rate",
                                    AllSelectionList.size(), 0,
                                    AllSelectionList.size());

  PThetaLep_outlier_low =
      SelectionHistsFromTOML<TH2D>("PThetaLep", plots_config);
  PThetaLep_outlier_low->SetName("PThetaLep_outlier_low");
  PThetaLep_outlier_low->SetZAxisTitle("Count with w <= 0.1");
  PThetaLep_outlier_high =
      SelectionHistsFromTOML<TH2D>("PThetaLep", plots_config);
  PThetaLep_outlier_high->SetName("PThetaLep_outlier_high");
  PThetaLep_outlier_high->SetZAxisTitle("Count with w >= 10");
  EnuQ2_outlier_low = SelectionHistsFromTOML<TH2D>("EnuQ2", plots_config);
  EnuQ2_outlier_low->SetName("EnuQ2_outlier_low");
  EnuQ2_outlier_low->SetZAxisTitle("Count with w <= 0.1");
  EnuQ2_outlier_high = SelectionHistsFromTOML<TH2D>("EnuQ2", plots_config);
  EnuQ2_outlier_high->SetName("EnuQ2_outlier_high");
  EnuQ2_outlier_high->SetZAxisTitle("Count with w >= 10");
  q0q3_high_outlier_low =
      SelectionHistsFromTOML<TH2D>("q0q3_high", plots_config);
  q0q3_high_outlier_low->SetName("q0q3_high_outlier_low");
  q0q3_high_outlier_low->SetZAxisTitle("Count with w <= 0.1");
  q0q3_high_outlier_high =
      SelectionHistsFromTOML<TH2D>("q0q3_high", plots_config);
  q0q3_high_outlier_high->SetName("q0q3_high_outlier_high");
  q0q3_high_outlier_high->SetZAxisTitle("Count with w >= 10");

  T2KNOvATruthTreeReader rdr(ttrdr);

  size_t nents = ttrdr.GetEntries(true);
  size_t ent_it = 0;
  size_t shout_it = nents / 25;

  double EnuCut = toml::find<double>(plots_config, "EnuCut");

  while (ttrdr.Next()) {

    int mode = bymode ? rdr.Mode() : 0;

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

    int primary_selection = rdr.GetPrimarySelection();
    double rw_w = 1;
    if (weightconfig == t2knova::kT2KND_to_NOvA) {
      rw_w *= t2knova::GetFakeDataWeight_ND280ToNOvA(
          rdr.PDGNu(), rdr.tgta(), rdr.Enu_true(), rdr.PLep(), rdr.AngLep_deg(),
          primary_selection);
    } else if (weightconfig == t2knova::kT2KND_to_T2KNonQE) {
      rw_w *= t2knova::GetFakeDataWeight_ND280ToT2KNonQE(
          rdr.PDGNu(), rdr.tgta(), rdr.Enu_true(), rdr.PLep(), rdr.AngLep_deg(),
          primary_selection);
    } else if (weightconfig == t2knova::kT2KND_to_T2KMnv1Pi) {
      rw_w *= t2knova::GetFakeDataWeight_ND280ToT2KMnv1Pi(
          rdr.PDGNu(), rdr.tgta(), rdr.Enu_true(), rdr.PLep(), rdr.AngLep_deg(),
          primary_selection);
    } else if (weightconfig == t2knova::kNOvA_to_T2KND_ptlep) {
      rw_w *= t2knova::GetFakeDataWeight_NOvAToT2KND_PtLep(
          rdr.PDGNu(), rdr.PDGLep(), rdr.tgta(), rdr.Enu_true(),
          (rdr.PLep()) * sqrt(1 - pow(rdr.CosLep(), 2)), rdr.Eav_NOvA(),
          primary_selection);
    } else if (weightconfig == t2knova::kNOvA_to_T2KPre_ptlep) {
      rw_w *= t2knova::GetFakeDataWeight_NOvAToT2KPre_PtLep(
          rdr.PDGNu(), rdr.PDGLep(), rdr.tgta(), rdr.Enu_true(),
          (rdr.PLep()) * sqrt(1 - pow(rdr.CosLep(), 2)), rdr.Eav_NOvA(),
          primary_selection);
    } else if (weightconfig == t2knova::kNOvA_to_T2KMnv1Pi_ptlep) {
      rw_w *= t2knova::GetFakeDataWeight_NOvAToT2KMnv1Pi_PtLep(
          rdr.PDGNu(), rdr.PDGLep(), rdr.tgta(), rdr.Enu_true(),
          (rdr.PLep()) * sqrt(1 - pow(rdr.CosLep(), 2)), rdr.Eav_NOvA(),
          primary_selection);
    } else if (weightconfig == t2knova::kNOvA_to_T2KNonQE_ptlep) {
      rw_w *= t2knova::GetFakeDataWeight_NOvAToT2KNonQE_PtLep(
          rdr.PDGNu(), rdr.PDGLep(), rdr.tgta(), rdr.Enu_true(),
          (rdr.PLep()) * sqrt(1 - pow(rdr.CosLep(), 2)), rdr.Eav_NOvA(),
          primary_selection);
    }
    w *= rw_w;

    Enu->Fill(w, sels, mode, rdr.Enu_true());
    ERecQE->Fill(w, sels, mode, rdr.GetTrueEnuQE());
    EnuERecQE->Fill(w, sels, mode, rdr.Enu_true(), rdr.GetTrueEnuQE());

    float ERecQEBias = (rdr.GetTrueEnuQE() / rdr.Enu_true()) - 1;
    EnuERecQEBias->Fill(w, sels, mode, rdr.Enu_true(), ERecQEBias);
    float ERecAvBias =
        ((rdr.Eav_NOvA() + rdr.FSLepP4().E()) / rdr.Enu_true()) - 1;
    EnuERecAvBias->Fill(w, sels, mode, rdr.Enu_true(), ERecAvBias);

    float q0_QE = rdr.GetTrueEnuQE() - rdr.FSLepP4().E();
    float EAvBias = (rdr.Eav_NOvA() / q0_QE) - 1;
    q0_QEEAvBias->Fill(w, sels, mode,  EAvBias, q0_QE);

    PLep->Fill(w, sels, mode, rdr.PLep());
    ThetaLep->Fill(w, sels, mode, rdr.AngLep_deg());
    CosThetaLep->Fill(w, sels, mode, rdr.CosLep());
    PThetaLep->Fill(w, sels, mode, rdr.PLep(), rdr.AngLep_deg());
    EAvHad->Fill(w, sels, mode, rdr.Eav_NOvA());
    PtLep->Fill(w, sels, mode, rdr.PLep() * sqrt(1 - pow(rdr.CosLep(), 2)));
    Q2->Fill(w, sels, mode, rdr.Q2());

    q0->Fill(w, sels, mode, rdr.q0());
    q3->Fill(w, sels, mode, rdr.q3());
    float yrecf = rdr.Eav_NOvA() / (rdr.Eav_NOvA() + rdr.FSLepP4().E());
    yrec->Fill(w, sels, mode, yrecf);
    Enuyrec->Fill(w, sels, mode, rdr.Enu_true(), yrecf);

    EnuQ2->Fill(w, sels, mode, rdr.Enu_true(), rdr.Q2());
    q0q3_low->Fill(w, sels, mode, rdr.q3(), rdr.q0());
    q0q3_high->Fill(w, sels, mode, rdr.q3(), rdr.q0());
    EnuEAvHad->Fill(w, sels, mode, rdr.Enu_true(), rdr.Eav_NOvA());
    hmfscpip->Fill(w, sels, mode, rdr.hmfscpip());
    hmfspi0p->Fill(w, sels, mode, rdr.hmfspi0p());
    ncpi->Fill(w, sels, mode, rdr.ncpi());
    npi0->Fill(w, sels, mode, rdr.npi0());
    hmfsprotonp->Fill(w, sels, mode, rdr.hmfsprotonp());
    hmfsneutronp->Fill(w, sels, mode, rdr.hmfsneutronp());
    nproton->Fill(w, sels, mode, rdr.nproton());
    nneutron->Fill(w, sels, mode, rdr.nneutron());
    EGamma->Fill(w, sels, mode, rdr.EGamma());
    EGamma_DeExcite->Fill(w, sels, mode, rdr.EGamma());

    for (auto sel : sels) {
      XSecs->Fill(w, mode, sel);
    }

    if (rw_w <= 0.1) {
      PThetaLep_outlier_low->Fill(1, sels, mode, rdr.PLep(), rdr.AngLep_deg());
      EnuQ2_outlier_low->Fill(1, sels, mode, rdr.Enu_true(), rdr.Q2());
      q0q3_high_outlier_low->Fill(1, sels, mode, rdr.q3(), rdr.q0());
    } else if (rw_w >= 10) {
      PThetaLep_outlier_high->Fill(1, sels, mode, rdr.PLep(), rdr.AngLep_deg());
      EnuQ2_outlier_high->Fill(1, sels, mode, rdr.Enu_true(), rdr.Q2());
      q0q3_high_outlier_high->Fill(1, sels, mode, rdr.q3(), rdr.q0());
    }

    ent_it++;
  }
}

std::string input_file = "", output_file = "", output_dir = "";
std::string hist_config_file = "";
std::string FDSInputs = "";

t2knova::reweightconfig wconfig = t2knova::kNoWeight;
int tgta_select = 0;

bool FromGenerated = true;

std::map<reweightconfig, std::string> inputhistnames;

std::map<reweightconfig, std::string> inputhistnames_FromGenerated = {
    {kT2KND_to_NOvA, "Generated_to_2020/EnuPLepThetaLep"},
    {kT2KND_to_T2KNonQE, "Generated_to_NonQE/EnuPLepThetaLep"},
    {kT2KND_to_T2KMnv1Pi, "Generated_to_Mnv1Pi/EnuPLepThetaLep"},
    {kNOvA_to_T2KND_ptlep, "Generated_to_BANFF_POST/EnuPtLepEAvHad"},
    {kNOvA_to_T2KPre_ptlep, "Generated_to_BANFF_PRE/EnuPtLepEAvHad"},
    {kNOvA_to_T2KMnv1Pi_ptlep, "Generated_to_Mnv1Pi/EnuPtLepEAvHad"},
    {kNOvA_to_T2KNonQE_ptlep, "Generated_to_NonQE/EnuPtLepEAvHad"},
};

std::map<reweightconfig, std::string> inputhistnames_FromBANFF_POST = {
    {kT2KND_to_NOvA, "BANFF_POST_to_2020/EnuPLepThetaLep"},
    {kT2KND_to_T2KNonQE, "BANFF_POST_to_NonQE/EnuPLepThetaLep"},
    {kT2KND_to_T2KMnv1Pi, "BANFF_POST_to_Mnv1Pi/EnuPLepThetaLep"},
    {kNOvA_to_T2KND_ptlep, "2020_to_BANFF_POST/EnuPtLepEAvHad"},
    {kNOvA_to_T2KPre_ptlep, "2020_to_BANFF_PRE/EnuPtLepEAvHad"},
    {kNOvA_to_T2KMnv1Pi_ptlep, "2020_to_Mnv1Pi/EnuPtLepEAvHad"},
    {kNOvA_to_T2KNonQE_ptlep, "2020_to_NonQE/EnuPtLepEAvHad"},
};

std::map<reweightconfig, std::string> inputhistnames_FromBANFF_PRE = {
    {kT2KND_to_NOvA, "BANFF_PRE_to_2020/EnuPLepThetaLep"},
    {kT2KND_to_T2KNonQE, "BANFF_PRE_to_NonQE/EnuPLepThetaLep"},
    {kT2KND_to_T2KMnv1Pi, "BANFF_PRE_to_Mnv1Pi/EnuPLepThetaLep"},
    {kNOvA_to_T2KND_ptlep, "2020_to_BANFF_POST/EnuPtLepEAvHad"},
    {kNOvA_to_T2KPre_ptlep, "2020_to_BANFF_PRE/EnuPtLepEAvHad"},
    {kNOvA_to_T2KMnv1Pi_ptlep, "2020_to_Mnv1Pi/EnuPtLepEAvHad"},
    {kNOvA_to_T2KNonQE_ptlep, "2020_to_NonQE/EnuPtLepEAvHad"},
};

void SayUsage(char const *argv[]) {
  std::cout << "[USAGE]: " << argv[0]
            << "\n"
               "\n\t-i <inp.root>            : Input file"
               "\n\t-F <FDSInputs.root>      : FDS Weighting file"
               "\n\t-H <config.toml>         : RWHist config file"
               "\n\t-o <out.root>            : Output file"
               "\n\t-M                       : Separate by Mode"
               "\n\t--oscillate              : Apply oscillation weights"
               "\n\t--From <BANFFPre|BANFFPost|Generated> : Weight as if "
               "Tuned/Not Tuned"
               "\n\t-a <[C|H|O|any]>         : Target descriptor"
               "\n\t-W <Config>              : ReWeight Config"
               "\n\t         Configs:"
               "\n\t              * T2KND_to_NOvA"
               "\n\t              * NOvA_to_T2KND_ptlep"
               "\n\t              * NOvA_to_T2KPre_ptlep"
               "\n\t              * NOvA_to_T2KMnv1Pi_ptlep"
               "\n\t              * NOvA_to_T2KNonQE_ptlep"
               "\n\t-T <Generated|NDTuned|Mnv1Pi|NonQE> : Tune to Apply."
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
    } else if (std::string(argv[opt]) == "--From") {
      std::string arg = std::string(argv[++opt]);
      inputhistnames = inputhistnames_FromGenerated;
      if (arg == "BANFFPre") {
        inputhistnames = inputhistnames_FromBANFF_PRE;
        FromGenerated = false;
      } else if (arg == "BANFFPost") {
        inputhistnames = inputhistnames_FromBANFF_POST;
        FromGenerated = false;
      } else if (arg != "Generated") {
        std::cout << "[ERROR]: Invalid option passed to --From, should be "
                     "BANFFPre, BANFFPost, or Generated."
                  << std::endl;
        abort();
      }

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
    } else if (std::string(argv[opt]) == "-T") {
      std::string arg = std::string(argv[++opt]);
      if (arg == "Generated") {
        FDSSet = kGenerated;
        std::cout << "Applying Generated Tune. " << std::endl;
      } else if (arg == "NDTuned") {
        FDSSet = kNDTuned;
        std::cout << "Applying NDTuned Tune. " << std::endl;
      } else if (arg == "Mnv1Pi") {
        FDSSet = kMnv1Pi;
        std::cout << "Applying Mnv1Pi Tune. " << std::endl;
      } else if (arg == "NonQE") {
        FDSSet = kNonQE;
        std::cout << "Applying NonQE Tune. " << std::endl;
      } else {
        std::cout << "Invalid FDS selector passed: " << argv[4]
                  << ". Should be None/Mnv1Pi/NonQE." << std::endl;
        abort();
      }
    } else if (std::string(argv[opt]) == "-W") {
      std::string arg = std::string(argv[++opt]);
      if (arg == "T2KND_to_NOvA") {
        wconfig = t2knova::kT2KND_to_NOvA;
      } else if (arg == "kT2KND_to_T2KNonQE") {
        wconfig = t2knova::kT2KND_to_T2KNonQE;
      } else if (arg == "kT2KND_to_T2KMnv1Pi") {
        wconfig = t2knova::kT2KND_to_T2KMnv1Pi;
      } else if (arg == "NOvA_to_T2KND_ptlep") {
        wconfig = t2knova::kNOvA_to_T2KND_ptlep;
      } else if (arg == "NOvA_to_T2KPre_ptlep") {
        wconfig = t2knova::kNOvA_to_T2KPre_ptlep;
      } else if (arg == "NOvA_to_T2KMnv1Pi_ptlep") {
        wconfig = t2knova::kNOvA_to_T2KMnv1Pi_ptlep;
      } else if (arg == "NOvA_to_T2KNonQE_ptlep") {
        wconfig = t2knova::kNOvA_to_T2KNonQE_ptlep;
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

  gStyle->SetOptStat(false);

  handleOpts(argc, argv);

  if ((wconfig != t2knova::kNoWeight) && !FDSInputs.length()) {
    std::cout << "[ERROR]: Asked to perform weighting, but no inputs passed."
              << std::endl;
    return 1;
  }

  if (FDSInputs.length()) {
    t2knova::LoadHists(FDSInputs, inputhistnames);
  }

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

  XSecs->Apply([](TH1D &h) {
    for (int i = 0; i < SelectionList.size(); ++i) {
      h.GetXaxis()->SetBinLabel(i + 1, SelectionList[i].c_str());
    }
  });

  XSecs->Write(dout, false);

  Enu->Write(dout, true);
  ERecQE->Write(dout, true);
  PLep->Write(dout, true);
  ThetaLep->Write(dout, true);
  CosThetaLep->Write(dout, true);
  PThetaLep->Write(dout, true);
  PThetaLep_outlier_low->Write(dout, true);
  PThetaLep_outlier_high->Write(dout, true);

  EAvHad->Write(dout, true);
  PtLep->Write(dout, true);
  Q2->Write(dout, true);
  EnuQ2->Write(dout, true);
  EnuQ2_outlier_low->Write(dout, true);
  EnuQ2_outlier_high->Write(dout, true);

  EnuERecQE->Write(dout, true);
  EnuERecQEBias->Write(dout, true);
  EnuERecAvBias->Write(dout, true);
  q0_QEEAvBias->Write(dout, true);
  q0q3_low->Write(dout, true);
  q0q3_high->Write(dout, true);
  q0q3_high_outlier_low->Write(dout, true);
  q0q3_high_outlier_high->Write(dout, true);
  EnuEAvHad->Write(dout, true);

  q0->Write(dout, true);
  q3->Write(dout, true);
  yrec->Write(dout, true);
  Enuyrec->Write(dout, true);

  hmfscpip->Write(dout, true);
  hmfspi0p->Write(dout, true);
  ncpi->Write(dout, true);
  npi0->Write(dout, true);

  hmfsprotonp->Write(dout, true);
  hmfsneutronp->Write(dout, true);
  nproton->Write(dout, true);
  nneutron->Write(dout, true);

  EGamma->Write(dout, true);
  EGamma_DeExcite->Write(dout, true);

  fout.Close();
}
