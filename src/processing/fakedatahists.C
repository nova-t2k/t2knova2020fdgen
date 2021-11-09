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

std::string FixProjName(std::string, std::string const &);

using namespace t2knova;

TrueChannelHist<TH1F> *totxsecs;
TrueChannelHist<TH1F> *totxsecs_untuned;

SelectionHists<TH1F> *Enu;
SelectionHists<TH1F> *Enu_untuned;
SelectionHists<TH1F> *Q2;
SelectionHists<TH1F> *Q2_untuned;
SelectionHists<TH1F> *EGamma;
SelectionHists<TH1F> *EGamma_untuned;
SelectionHists<TH3F> *EnuPLepThetaLep;
SelectionHists<TH3F> *EnuPLepThetaLep_untuned;
SelectionHists<TH3F> *EnuPLepEAvHad;
SelectionHists<TH3F> *EnuQ2EAvHad;
SelectionHists<TH3F> *EnuPtLepEAvHad;

void Fill(TTreeReader &ttrdr, toml::value const &plots_config, bool ist2k,
          int tgta_select = 0) {
  Enu = SelectionHistsFromTOML<TH1F>("Enu", toml::find(plots_config, "Enu"));
  Enu_untuned = SelectionHistsFromTOML<TH1F>("Enu_untuned",
                                             toml::find(plots_config, "Enu"));

  Q2 = SelectionHistsFromTOML<TH1F>("Q2", toml::find(plots_config, "Q2"));
  Q2_untuned = SelectionHistsFromTOML<TH1F>("Q2_untuned",
                                            toml::find(plots_config, "Q2"));

  EGamma = SelectionHistsFromTOML<TH1F>("EGamma",
                                        toml::find(plots_config, "EGamma"));
  EGamma_untuned = SelectionHistsFromTOML<TH1F>(
      "EGamma_untuned", toml::find(plots_config, "EGamma"));

  if (ist2k) {
    EnuPLepThetaLep = SelectionHistsFromTOML<TH3F>(
        "EnuPLepThetaLep", toml::find(plots_config, "EnuPLepThetaLep"));
    EnuPLepThetaLep_untuned = SelectionHistsFromTOML<TH3F>(
        "EnuPLepThetaLep_untuned", toml::find(plots_config, "EnuPLepThetaLep"));
  } else {
    EnuPLepEAvHad = SelectionHistsFromTOML<TH3F>(
        "EnuPLepEAvHad", toml::find(plots_config, "EnuPLepEAvHad"));
    EnuQ2EAvHad = SelectionHistsFromTOML<TH3F>(
        "EnuQ2EAvHad", toml::find(plots_config, "EnuQ2EAvHad"));
    EnuPtLepEAvHad = SelectionHistsFromTOML<TH3F>(
        "EnuPtLepEAvHad", toml::find(plots_config, "EnuPtLepEAvHad"));
  }

  totxsecs =
      new TrueChannelHist<TH1F>("totxsecs", ";;#sigma^{#int#Phi} cm^{2}",
                                SelectionList.size(), 0, SelectionList.size());

  totxsecs_untuned = new TrueChannelHist<TH1F>(
      "totxsecs_untuned", ";;#sigma^{#int#Phi} cm^{2}", SelectionList.size(), 0,
      SelectionList.size());

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

    bool iscc = rdr.PDGLep() % 2;

    for (int s : sels) {
      totxsecs->Fill(w, rdr.Mode(), s);
      totxsecs_untuned->Fill(rdr.fScaleFactor(), rdr.Mode(), s);
    }

    Enu->Fill(w, sels, rdr.Mode(), rdr.Enu_true());
    Enu_untuned->Fill(rdr.fScaleFactor(), sels, rdr.Mode(), rdr.Enu_true());
    Q2->Fill(w, sels, rdr.Mode(), rdr.Q2());
    Q2_untuned->Fill(rdr.fScaleFactor(), sels, rdr.Mode(), rdr.Q2());
    EGamma->Fill(w, sels, rdr.Mode(), rdr.EGamma());

    EGamma_untuned->Fill(rdr.fScaleFactor(), sels, rdr.Mode(), rdr.EGamma());
    if (ist2k) {
      EnuPLepThetaLep->Fill(w, sels, rdr.Mode(), rdr.Enu_true(), rdr.PLep(),
                            rdr.AngLep_deg());
      EnuPLepThetaLep_untuned->Fill(rdr.fScaleFactor(), sels, rdr.Mode(),
                                    rdr.Enu_true(), rdr.PLep(),
                                    rdr.AngLep_deg());
    } else {
      EnuPLepEAvHad->Fill(w, sels, rdr.Mode(), rdr.Enu_true(), rdr.PLep(),
                          rdr.EavAlt());
      EnuQ2EAvHad->Fill(w, sels, rdr.Mode(), rdr.Enu_true(), rdr.Q2(),
                        rdr.EavAlt());
      EnuPtLepEAvHad->Fill(w, sels, rdr.Mode(), rdr.Enu_true(),
                           rdr.PLep() * sqrt(1 - pow(rdr.CosLep(), 2)),
                           rdr.EavAlt());
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

  auto axis_labeler = [](TH1F &h) {
    for (int i = 0; i < SelectionList.size(); ++i) {
      h.GetXaxis()->SetBinLabel(i + 1, SelectionList[i].c_str());
    }
  };

  totxsecs->Apply(axis_labeler);
  totxsecs->Write(dout);

  totxsecs_untuned->Apply(axis_labeler);
  totxsecs_untuned->Write(dout);

  Enu->Write(dout, true);
  Enu_untuned->Write(dout, true);

  Q2->Write(dout, true);
  Q2_untuned->Write(dout, true);

  EGamma->Write(dout, true);
  EGamma_untuned->Write(dout, true);

  auto ProjYZ = [=](TH3F const &h) -> TH2F {
    std::unique_ptr<TH2F> h2 = dynamic_cast_uptr<TH2F>(Project3D(h, "yz"));
    h2->SetDirectory(nullptr);
    h2->SetName(FixProjName(h.GetName(), "yz").c_str());
    return TH2F(*h2.get());
  };
  auto ProjYX = [=](TH3F const &h) -> TH2F {
    std::unique_ptr<TH2F> h2 = dynamic_cast_uptr<TH2F>(Project3D(h, "yx"));
    h2->SetDirectory(nullptr);
    h2->SetName(FixProjName(h.GetName(), "yx").c_str());
    return TH2F(*h2.get());
  };
  auto ProjXZ = [=](TH3F const &h) -> TH2F {
    std::unique_ptr<TH2F> h2 = dynamic_cast_uptr<TH2F>(Project3D(h, "xz"));
    h2->SetDirectory(nullptr);
    h2->SetName(FixProjName(h.GetName(), "xz").c_str());
    return TH2F(*h2.get());
  };

  auto ProjX = [=](TH3F const &h) -> TH1F {
    std::unique_ptr<TH1F> h1 = dynamic_cast_uptr<TH1F>(Project3D(h, "x"));
    h1->SetDirectory(nullptr);
    h1->SetName(FixProjName(h.GetName(), "x").c_str());
    return TH1F(*h1.get());
  };

  auto ProjY = [=](TH3F const &h) -> TH1F {
    std::unique_ptr<TH1F> h1 = dynamic_cast_uptr<TH1F>(Project3D(h, "y"));
    h1->SetDirectory(nullptr);
    h1->SetName(FixProjName(h.GetName(), "y").c_str());
    return TH1F(*h1.get());
  };

  auto ProjZ = [=](TH3F const &h) -> TH1F {
    std::unique_ptr<TH1F> h1 = dynamic_cast_uptr<TH1F>(Project3D(h, "z"));
    h1->SetDirectory(nullptr);
    h1->SetName(FixProjName(h.GetName(), "z").c_str());
    return TH1F(*h1.get());
  };

  if (ist2k) {
    EnuPLepThetaLep->Write(dout, true);
    auto PLepThetaLep = EnuPLepThetaLep->Transform<TH2F>(ProjYZ);
    PLepThetaLep.Write(dout, true);
    auto EnuPLep = EnuPLepThetaLep->Transform<TH2F>(ProjYX);
    EnuPLep.Write(dout, true);
    auto EnuThetaLep = EnuPLepThetaLep->Transform<TH2F>(ProjXZ);
    EnuThetaLep.Write(dout, true);

    auto EnuProj = EnuPLepThetaLep->Transform<TH1F>(ProjX);
    EnuProj.Write(dout, true);
    auto PLep = EnuPLepThetaLep->Transform<TH1F>(ProjY);
    PLep.Write(dout, true);
    auto ThetaLep = EnuPLepThetaLep->Transform<TH1F>(ProjZ);
    ThetaLep.Write(dout, true);

    EnuPLepThetaLep_untuned->Write(dout, true);
    auto PLepThetaLep_untuned =
        EnuPLepThetaLep_untuned->Transform<TH2F>(ProjYZ);
    PLepThetaLep_untuned.Write(dout, true);
    auto EnuPLep_untuned = EnuPLepThetaLep_untuned->Transform<TH2F>(ProjYX);
    EnuPLep_untuned.Write(dout, true);
    auto EnuThetaLep_untuned = EnuPLepThetaLep_untuned->Transform<TH2F>(ProjXZ);
    EnuThetaLep_untuned.Write(dout, true);

    auto EnuProj_untuned = EnuPLepThetaLep_untuned->Transform<TH1F>(ProjX);
    EnuProj_untuned.Write(dout, true);
    auto PLep_untuned = EnuPLepThetaLep_untuned->Transform<TH1F>(ProjY);
    PLep_untuned.Write(dout, true);
    auto ThetaLep_untuned = EnuPLepThetaLep_untuned->Transform<TH1F>(ProjZ);
    ThetaLep_untuned.Write(dout, true);
  } else {
    EnuPLepEAvHad->Write(dout, true);
    auto PLepEAvHad = EnuPLepEAvHad->Transform<TH2F>(ProjYZ);
    PLepEAvHad.Write(dout, true);
    auto EnuPLep = EnuPLepEAvHad->Transform<TH2F>(ProjYX);
    EnuPLep.Write(dout, true);
    auto EnuEAvHad = EnuPLepEAvHad->Transform<TH2F>(ProjXZ);
    EnuEAvHad.Write(dout, true);

    EnuQ2EAvHad->Write(dout, true);
    auto Q2EAvHad = EnuQ2EAvHad->Transform<TH2F>(ProjYZ);
    Q2EAvHad.Write(dout, true);
    auto EnuQ2 = EnuQ2EAvHad->Transform<TH2F>(ProjYX);
    EnuQ2.Write(dout, true);

    EnuPtLepEAvHad->Write(dout, true);
    auto PtLepEAvHad = EnuPtLepEAvHad->Transform<TH2F>(ProjYZ);
    PtLepEAvHad.Write(dout, true);
    auto EnuPtLep = EnuPtLepEAvHad->Transform<TH2F>(ProjYX);
    EnuPtLep.Write(dout, true);
  }

  fout.Close();
}

std::string FixProjName(std::string name, std::string const &proj) {
  std::string ogname = name;
  if (name.find("EnuPLepThetaLep") != std::string::npos) {
    if (proj == "yz") {
      name.replace(0, std::string("EnuPLepThetaLep").size(),
                   "PLepThetaLep_prj0");
    } else if (proj == "yx") {
      name.replace(0, std::string("EnuPLepThetaLep").size(), "PLepEnu_prj0");
    } else if (proj == "xz") {
      name.replace(0, std::string("EnuPLepThetaLep").size(),
                   "EnuThetaLep_prj0");
    } else if (proj == "x") {
      name.replace(0, std::string("EnuPLepThetaLep").size(), "Enu_prj0");
    } else if (proj == "y") {
      name.replace(0, std::string("EnuPLepThetaLep").size(), "PLep_prj0");
    } else if (proj == "z") {
      name.replace(0, std::string("EnuPLepThetaLep").size(), "ThetaLep_prj0");
    }
  } else if (name.find("EnuPLepEAvHad") != std::string::npos) {
    if (proj == "yz") {
      name.replace(0, std::string("EnuPLepEAvHad").size(), "PLepEAvHad_prj1");
    } else if (proj == "yx") {
      name.replace(0, std::string("EnuPLepEAvHad").size(), "PLepEnu_prj1");
    } else if (proj == "xz") {
      name.replace(0, std::string("EnuPLepEAvHad").size(), "EnuEAvHad_prj1");
    } else if (proj == "x") {
      name.replace(0, std::string("EnuPLepEAvHad").size(), "Enu_prj1");
    } else if (proj == "y") {
      name.replace(0, std::string("EnuPLepEAvHad").size(), "PLep_prj1");
    } else if (proj == "z") {
      name.replace(0, std::string("EnuPLepEAvHad").size(), "EAvHad_prj1");
    }
  } else if (name.find("EnuQ2EAvHad") != std::string::npos) {
    if (proj == "yz") {
      name.replace(0, std::string("EnuQ2EAvHad").size(), "Q2EAvHad_prj2");
    } else if (proj == "yx") {
      name.replace(0, std::string("EnuQ2EAvHad").size(), "Q2Enu_prj2");
    } else if (proj == "xz") {
      name.replace(0, std::string("EnuQ2EAvHad").size(), "EnuEAvHad_prj2");
    } else if (proj == "x") {
      name.replace(0, std::string("EnuQ2EAvHad").size(), "Enu_prj2");
    } else if (proj == "y") {
      name.replace(0, std::string("EnuQ2EAvHad").size(), "Q2_prj2");
    } else if (proj == "z") {
      name.replace(0, std::string("EnuQ2EAvHad").size(), "EAvHad_prj2");
    }
  } else if (name.find("EnuPtLepEAvHad") != std::string::npos) {
    if (proj == "yz") {
      name.replace(0, std::string("EnuPtLepEAvHad").size(), "PtLepEAvHad_prj3");
    } else if (proj == "yx") {
      name.replace(0, std::string("EnuPtLepEAvHad").size(), "PtLepEnu_prj3");
    } else if (proj == "xz") {
      name.replace(0, std::string("EnuPtLepEAvHad").size(), "EnuEAvHad_prj3");
    } else if (proj == "x") {
      name.replace(0, std::string("EnuPtLepEAvHad").size(), "Enu_prj3");
    } else if (proj == "y") {
      name.replace(0, std::string("EnuPtLepEAvHad").size(), "PtLep_prj3");
    } else if (proj == "z") {
      name.replace(0, std::string("EnuPtLepEAvHad").size(), "EAvHad_prj3");
    }
  }
  return name;
}