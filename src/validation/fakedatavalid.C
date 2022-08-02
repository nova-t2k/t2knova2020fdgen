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

#include "OscillationHelper.hxx"

#include <cmath>

using namespace t2knova;

bool bymode = false;
bool doosc = false;

enum FDS { kGenerated, kNDTuned, kMnv1Pi, kNonQE };
FDS FDSSet = kGenerated;

TrueChannelHist<TH1F> *XSecs;
SelectionHists<TH1F> *EvWeights;

SelectionHists<TH1F> *Enu;
SelectionHists<TH1F> *ERecQE;
SelectionHists<TH1F> *PLep;
// SelectionHists<TH1F> *ThetaLep;
SelectionHists<TH1F> *CosThetaLep;
SelectionHists<TH2F> *PThetaLep;
// SelectionHists<TH2F> *PThetaLep_outlier_low;
// SelectionHists<TH2F> *PThetaLep_outlier_high;
SelectionHists<TH1F> *EAvHad;
SelectionHists<TH1F> *PtLep;
SelectionHists<TH1F> *Q2;
SelectionHists<TH1F> *q0;
SelectionHists<TH1F> *yrec;
SelectionHists<TH2F> *Enuyrec;
SelectionHists<TH2F> *EnuQ2;
// SelectionHists<TH2F> *EnuQ2_outlier_low;
// SelectionHists<TH2F> *EnuQ2_outlier_high;
SelectionHists<TH2F> *EnuERecQE;
SelectionHists<TH2F> *EnuERecQEBias;
SelectionHists<TH2F> *EnuERecAvBias;
// SelectionHists<TH2F> *q0q3_high_outlier_low;
// SelectionHists<TH2F> *q0q3_high_outlier_high;
SelectionHists<TH2F> *q0q3_low;
SelectionHists<TH2F> *q0q3_high;
SelectionHists<TH2F> *EnuEAvHad;
SelectionHists<TH1F> *hmfscpip;
SelectionHists<TH1F> *hmfspi0p;
SelectionHists<TH1F> *ncpi;
SelectionHists<TH1F> *npi0;
SelectionHists<TH1F> *hmfsprotonp;
SelectionHists<TH1F> *hmfsneutronp;
SelectionHists<TH1F> *nproton;
SelectionHists<TH1F> *nneutron;
// SelectionHists<TH1F> *EGamma;
// SelectionHists<TH1F> *EGamma_DeExcite;

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

  Enu = SelectionHistsFromTOML<TH1F>("Enu", plots_config);

  ERecQE = SelectionHistsFromTOML<TH1F>("ERecQE", plots_config);
  PLep = SelectionHistsFromTOML<TH1F>("PLep", plots_config);
  Q2 = SelectionHistsFromTOML<TH1F>("Q2", plots_config);
  EnuQ2 = SelectionHistsFromTOML<TH2F>("EnuQ2", plots_config);
  EnuERecQE = SelectionHistsFromTOML<TH2F>("EnuERecQE", plots_config);
  EnuERecQEBias = SelectionHistsFromTOML<TH2F>("EnuERecQEBias", plots_config);
  EnuERecAvBias = SelectionHistsFromTOML<TH2F>("EnuERecAvBias", plots_config);

  q0q3_low = SelectionHistsFromTOML<TH2F>("q0q3_low", plots_config);
  q0q3_high = SelectionHistsFromTOML<TH2F>("q0q3_high", plots_config);
  EnuEAvHad = SelectionHistsFromTOML<TH2F>("EnuEAvHad", plots_config);

  q0 = SelectionHistsFromTOML<TH1F>("q0", plots_config);
  yrec = SelectionHistsFromTOML<TH1F>("yrec", plots_config);
  Enuyrec = SelectionHistsFromTOML<TH2F>("Enuyrec", plots_config);

  // ThetaLep = SelectionHistsFromTOML<TH1F>("ThetaLep", plots_config);
  CosThetaLep = SelectionHistsFromTOML<TH1F>("CosThetaLep", plots_config);

  PThetaLep = SelectionHistsFromTOML<TH2F>("PThetaLep", plots_config);

  EAvHad = SelectionHistsFromTOML<TH1F>("EAvHad", plots_config);

  PtLep = SelectionHistsFromTOML<TH1F>("PtLep", plots_config);

  hmfscpip = SelectionHistsFromTOML<TH1F>("hmfscpip", plots_config);
  hmfspi0p = SelectionHistsFromTOML<TH1F>("hmfspi0p", plots_config);
  ncpi = SelectionHistsFromTOML<TH1F>("ncpi", plots_config);
  npi0 = SelectionHistsFromTOML<TH1F>("npi0", plots_config);

  hmfsprotonp = SelectionHistsFromTOML<TH1F>("hmfsprotonp", plots_config);
  hmfsneutronp = SelectionHistsFromTOML<TH1F>("hmfsneutronp", plots_config);
  nproton = SelectionHistsFromTOML<TH1F>("nproton", plots_config);
  nneutron = SelectionHistsFromTOML<TH1F>("nneutron", plots_config);

  // EGamma = SelectionHistsFromTOML<TH1F>("EGamma", plots_config);
  // EGamma_DeExcite =
  //     SelectionHistsFromTOML<TH1F>("EGamma_DeExcite", plots_config);

  XSecs = new TrueChannelHist<TH1F>("SelectionXSecs", ";Selection;Rate",
                                    AllSelectionList.size(), 0,
                                    AllSelectionList.size());
  EvWeights = SelectionHistsFromTOML<TH1F>("EvWeights", plots_config);

  // PThetaLep_outlier_low =
  //     SelectionHistsFromTOML<TH2F>("PThetaLep", plots_config);
  // PThetaLep_outlier_low->SetName("PThetaLep_outlier_low");
  // PThetaLep_outlier_low->SetZAxisTitle("Count with w <= 0.1");
  // PThetaLep_outlier_high =
  //     SelectionHistsFromTOML<TH2F>("PThetaLep", plots_config);
  // PThetaLep_outlier_high->SetName("PThetaLep_outlier_high");
  // PThetaLep_outlier_high->SetZAxisTitle("Count with w >= 10");
  // EnuQ2_outlier_low = SelectionHistsFromTOML<TH2F>("EnuQ2", plots_config);
  // EnuQ2_outlier_low->SetName("EnuQ2_outlier_low");
  // EnuQ2_outlier_low->SetZAxisTitle("Count with w <= 0.1");
  // EnuQ2_outlier_high = SelectionHistsFromTOML<TH2F>("EnuQ2", plots_config);
  // EnuQ2_outlier_high->SetName("EnuQ2_outlier_high");
  // EnuQ2_outlier_high->SetZAxisTitle("Count with w >= 10");
  // q0q3_high_outlier_low =
  //     SelectionHistsFromTOML<TH2F>("q0q3_high", plots_config);
  // q0q3_high_outlier_low->SetName("q0q3_high_outlier_low");
  // q0q3_high_outlier_low->SetZAxisTitle("Count with w <= 0.1");
  // q0q3_high_outlier_high =
  //     SelectionHistsFromTOML<TH2F>("q0q3_high", plots_config);
  // q0q3_high_outlier_high->SetName("q0q3_high_outlier_high");
  // q0q3_high_outlier_high->SetZAxisTitle("Count with w >= 10");

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
          w *= GetnonQEWeight(rdr.PDGNu(), rdr.GetQ2QE());
        }
      }
    }


    if (doosc && (std::abs(rdr.Mode()) < 30)) { // Oscillate CC events if
                                                // enabled
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
    double rw_w = 1;
    if (weightconfig == t2knova::kT2KND_to_NOvA) {
      rw_w *= t2knova::GetFakeDataWeight_ND280ToNOvA(
          rdr.PDGNu(), rdr.PDGLep(), rdr.tgta(), rdr.Enu_true(), rdr.PLep(),
          rdr.AngLep_deg(), primary_selection);
    } else if (weightconfig == t2knova::kT2KND_to_T2KNonQE) {
      rw_w *= t2knova::GetFakeDataWeight_ND280ToT2KNonQE(
          rdr.PDGNu(), rdr.PDGLep(), rdr.tgta(), rdr.Enu_true(), rdr.PLep(),
          rdr.AngLep_deg(), primary_selection);
    } else if (weightconfig == t2knova::kT2KND_to_T2KMnv1Pi) {
      rw_w *= t2knova::GetFakeDataWeight_ND280ToT2KMnv1Pi(
          rdr.PDGNu(), rdr.PDGLep(), rdr.tgta(), rdr.Enu_true(), rdr.PLep(),
          rdr.AngLep_deg(), primary_selection);
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

    Enu->Fill(w, sels, rdr.Mode(), rdr.Enu_true());
    double erec = EnuQErec(rdr.FSLepP4().E(), rdr.PLep(), rdr.CosLep(), 0,
                           rdr.PDGNu() > 0);
    ERecQE->Fill(w, sels, rdr.Mode(), erec);
    EnuERecQE->Fill(w, sels, rdr.Mode(), rdr.Enu_true(), erec);

    float ERecQEBias = (erec / rdr.Enu_true()) - 1;
    EnuERecQEBias->Fill(w, sels, rdr.Mode(), rdr.Enu_true(), ERecQEBias);
    float ERecAvBias =
        ((rdr.Eav_NOvA() + rdr.FSLepP4().E()) / rdr.Enu_true()) - 1;
    EnuERecAvBias->Fill(w, sels, rdr.Mode(), rdr.Enu_true(), ERecAvBias);

    PLep->Fill(w, sels, rdr.Mode(), rdr.PLep());
    // ThetaLep->Fill(w, sels, rdr.Mode(), rdr.AngLep_deg());
    CosThetaLep->Fill(w, sels, rdr.Mode(), rdr.CosLep());
    PThetaLep->Fill(w, sels, rdr.Mode(), rdr.PLep(), rdr.AngLep_deg());
    EAvHad->Fill(w, sels, rdr.Mode(), rdr.Eav_NOvA());
    PtLep->Fill(w, sels, rdr.Mode(),
                rdr.PLep() * sqrt(1 - pow(rdr.CosLep(), 2)));
    Q2->Fill(w, sels, rdr.Mode(), rdr.Q2());

    q0->Fill(w, sels, rdr.Mode(), rdr.q0());
    float yrecf = rdr.Eav_NOvA() / (rdr.Eav_NOvA() + rdr.FSLepP4().E());
    yrec->Fill(w, sels, rdr.Mode(), yrecf);
    Enuyrec->Fill(w, sels, rdr.Mode(), rdr.Enu_true(), yrecf);

    EnuQ2->Fill(w, sels, rdr.Mode(), rdr.Enu_true(), rdr.Q2());
    q0q3_low->Fill(w, sels, rdr.Mode(), rdr.q3(), rdr.q0());
    q0q3_high->Fill(w, sels, rdr.Mode(), rdr.q3(), rdr.q0());
    EnuEAvHad->Fill(w, sels, rdr.Mode(), rdr.Enu_true(), rdr.Eav_NOvA());
    hmfscpip->Fill(w, sels, rdr.Mode(), rdr.hmfscpip());
    hmfspi0p->Fill(w, sels, rdr.Mode(), rdr.hmfspi0p());
    ncpi->Fill(w, sels, rdr.Mode(), rdr.ncpi());
    npi0->Fill(w, sels, rdr.Mode(), rdr.npi0());
    hmfsprotonp->Fill(w, sels, rdr.Mode(), rdr.hmfsprotonp());
    hmfsneutronp->Fill(w, sels, rdr.Mode(), rdr.hmfsneutronp());
    nproton->Fill(w, sels, rdr.Mode(), rdr.nproton());
    nneutron->Fill(w, sels, rdr.Mode(), rdr.nneutron());
    // EGamma->Fill(w, sels, rdr.Mode(), rdr.EGamma());
    // EGamma_DeExcite->Fill(w, sels, rdr.Mode(), rdr.EGamma());

    for (auto sel : sels) {
      XSecs->Fill(w, rdr.Mode(), sel);
    }

    EvWeights->Fill(1, sels, rdr.Mode(), w);

    // if (rw_w <= 0.1) {
    //   PThetaLep_outlier_low->Fill(1, sels, rdr.Mode(), rdr.PLep(),
    //                               rdr.AngLep_deg());
    //   EnuQ2_outlier_low->Fill(1, sels, rdr.Mode(), rdr.Enu_true(), rdr.Q2());
    //   q0q3_high_outlier_low->Fill(1, sels, rdr.Mode(), rdr.q3(), rdr.q0());
    // } else if (rw_w >= 10) {
    //   PThetaLep_outlier_high->Fill(1, sels, rdr.Mode(), rdr.PLep(),
    //                                rdr.AngLep_deg());
    //   EnuQ2_outlier_high->Fill(1, sels, rdr.Mode(), rdr.Enu_true(), rdr.Q2());
    //   q0q3_high_outlier_high->Fill(1, sels, rdr.Mode(), rdr.q3(), rdr.q0());
    // }

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
    {kT2KND_to_T2KMnv1Pi, "BANFF_PRE_to_NonQE/EnuPLepThetaLep"},
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
               "\n\t--No-Tune                : Don't apply tune weights"
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
    } else if (std::string(argv[opt]) == "-T") {
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
    }else if (std::string(argv[opt]) == "-W") {
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

  XSecs->Apply([](TH1F &h) {
    for (int i = 0; i < SelectionList.size(); ++i) {
      h.GetXaxis()->SetBinLabel(i + 1, SelectionList[i].c_str());
    }
  });

  XSecs->Write(dout, false);

  Enu->Write(dout, true);
  ERecQE->Write(dout, true);
  PLep->Write(dout, true);
  // ThetaLep->Write(dout, true);
  CosThetaLep->Write(dout, true);
  PThetaLep->Write(dout, true);
  // PThetaLep_outlier_low->Write(dout, true);
  // PThetaLep_outlier_high->Write(dout, true);

  EAvHad->Write(dout, true);
  PtLep->Write(dout, true);
  Q2->Write(dout, true);
  EnuQ2->Write(dout, true);
  // EnuQ2_outlier_low->Write(dout, true);
  // EnuQ2_outlier_high->Write(dout, true);

  EnuERecQE->Write(dout, true);
  EnuERecQEBias->Write(dout, true);
  EnuERecAvBias->Write(dout, true);
  // q0q3_low->Write(dout, true);
  // q0q3_high->Write(dout, true);
  // q0q3_high_outlier_low->Write(dout, true);
  // q0q3_high_outlier_high->Write(dout, true);
  EnuEAvHad->Write(dout, true);

  q0->Write(dout, true);
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

  // EGamma->Write(dout, true);
  // EGamma_DeExcite->Write(dout, true);

  EvWeights->Write(dout, true);

  fout.Close();
}
