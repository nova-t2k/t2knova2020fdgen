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

#include "hblob.h"

using namespace t2knova;

modeblob<TH1D> *totxsecs;
modeblob<TH1D> *totxsecs_untuned;

hblob<TH1D> *Enu;
hblob<TH1D> *Enu_untuned;
hblob<TH1D> *Q2;
hblob<TH1D> *Q2_untuned;
hblob<TH3D> *EnuPLepThetaLep;
hblob<TH3D> *EnuPLepThetaLep_untuned;
hblob<TH3D> *EnuPLepEAvHad;
hblob<TH3D> *EnuQ2EAvHad;
hblob<TH3D> *EnuPtLepEAvHad;

std::string FixProjName(std::string name, std::string const &proj) {

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
}

void Fill(TTreeReader &rdr, bool ist2k, int tgta_select = 0) {
  Enu = new hblob<TH1D>(
      "Enu",
      ";#it{E_{#nu}} (GeV);#it{y}; d#sigma/d#it{E_{#nu}} cm^{2} GeV^{-1}", 200,
      0, 10);
  Enu_untuned = new hblob<TH1D>(
      "Enu_untuned",
      ";#it{E_{#nu}} (GeV);#it{y}; d#sigma/d#it{E_{#nu}} cm^{2} GeV^{-1}", 200,
      0, 10);

  Q2 = new hblob<TH1D>(
      "Q2", ";#it{Q}^{2} (GeV^{2});#it{y}; d#sigma/d#it{Q}^{2} cm^{2} GeV^{-2}",
      200, 0, 4);
  Q2_untuned = new hblob<TH1D>(
      "Q2_untuned",
      ";#it{Q}^{2} (GeV^{2});#it{y}; d#sigma/d#it{Q}^{2} cm^{2} GeV^{-2}", 200,
      0, 4);

  if (ist2k) {
    EnuPLepThetaLep = new hblob<TH3D>(
        "EnuPLepThetaLep",
        ";#it{E_{#nu}} (GeV);#it{p}_{lep} (GeV/#it{c});#theta_{lep} (rad)"
        "GeV^{-2} #it{c} rad^{-1}",
        100, 0, 5, 50, 0, 5, 25, 0, M_PI);

    EnuPLepThetaLep_untuned = new hblob<TH3D>(
        "EnuPLepThetaLep_untuned",
        ";#it{E_{#nu}} (GeV);#it{p}_{lep} (GeV/#it{c});#theta_{lep} (rad)"
        "GeV^{-2} #it{c} rad^{-1}",
        100, 0, 5, 50, 0, 5, 25, 0, M_PI);
  } else {
    EnuPLepEAvHad = new hblob<TH3D>(
        "EnuPLepEAvHad",
        ";#it{E_{#nu}} (GeV);#it{p}_{lep} (GeV/#it{c}); #it{E}_{Had}^{Vis} "
        "(GeV)",
        100, 0, 5, 50, 0, 5, 40, 0, 4);
    EnuQ2EAvHad = new hblob<TH3D>(
        "EnuQ2EAvHad",
        ";#it{E_{#nu}} (GeV);Q^{2} (GeV^{2}/#it{c}^{2}); #it{E}_{Had}^{Vis} "
        "(GeV)",
        100, 0, 5, 50, 0, 5, 40, 0, 4);
    EnuPtLepEAvHad = new hblob<TH3D>(
        "EnuPtLepEAvHad",
        ";#it{E_{#nu}} (GeV);#it{p}^{T}_{lep} (GeV/#it{c}); #it{E}_{Had}^{Vis} "
        "(GeV)",
        100, 0, 5, 40, 0, 2, 40, 0, 4);
  }

  totxsecs =
      new modeblob<TH1D>("totxsecs", ";;#sigma^{#int#Phi} cm^{2}", 10, 0, 10);

  totxsecs_untuned = new modeblob<TH1D>(
      "totxsecs_untuned", ";;#sigma^{#int#Phi} cm^{2}", 10, 0, 10);

  TTreeReaderValue<int> PDGLep(rdr, "PDGLep");
  TTreeReaderValue<int> Mode(rdr, "Mode");

  TTreeReaderValue<float> CosLep(rdr, "CosLep");
  TTreeReaderValue<float> EavAlt(rdr, "EavAlt");
  TTreeReaderValue<float> Enu_true(rdr, "Enu_true");
  TTreeReaderValue<float> PLep_v(rdr, "PLep");
  TTreeReaderValue<float> Q2var(rdr, "Q2");
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

    bool iscc = (*PDGLep) % 2;

    int NFSCPi = *NFSpip + *NFSpim;
    int NFSOther_gammaextralep = *NFSOther + *NFSgamma;
    if ((*NFSlep) > 1) {
      NFSOther_gammaextralep += (*NFSlep) - 1;
    }

    FlagBlob fblob =
        GetFlagBlob(iscc, NFSCPi, *NFSpi0, NFSOther_gammaextralep, *Mode);

    if (fblob.flagCCINC) {
      totxsecs->Fill(w, *Mode, 0.0);
      totxsecs_untuned->Fill(*fScaleFactor, *Mode, 0.0);
      if (fblob.flagCC0pi) {
        totxsecs->Fill(w, *Mode, 1.0);
        totxsecs_untuned->Fill(*fScaleFactor, *Mode, 1.0);
      } else if (fblob.flagCC1cpi) {
        totxsecs->Fill(w, *Mode, 2.0);
        totxsecs_untuned->Fill(*fScaleFactor, *Mode, 2.0);
      } else if (fblob.flagCC1pi0) {
        totxsecs->Fill(w, *Mode, 3.0);
        totxsecs_untuned->Fill(*fScaleFactor, *Mode, 3.0);
      } else {
        totxsecs->Fill(w, *Mode, 4.0);
        totxsecs_untuned->Fill(*fScaleFactor, *Mode, 4.0);
      }
    } else if (fblob.flagNCINC) {
      totxsecs->Fill(w, *Mode, 5.0);
      totxsecs_untuned->Fill(*fScaleFactor, *Mode, 5.0);
      if (fblob.flagNC0pi) {
        totxsecs->Fill(w, *Mode, 6.0);
        totxsecs_untuned->Fill(*fScaleFactor, *Mode, 6.0);
      } else if (fblob.flagNC1cpi) {
        totxsecs->Fill(w, *Mode, 7.0);
        totxsecs_untuned->Fill(*fScaleFactor, *Mode, 7.0);
      } else if (fblob.flagNC1pi0) {
        totxsecs->Fill(w, *Mode, 8.0);
        totxsecs_untuned->Fill(*fScaleFactor, *Mode, 8.0);
      } else {
        totxsecs->Fill(w, *Mode, 9.0);
        totxsecs_untuned->Fill(*fScaleFactor, *Mode, 9.0);
      }
    }

    Enu->Fill(w, fblob, *Enu_true);
    Enu_untuned->Fill(*fScaleFactor, fblob, *Enu_true);
    Q2->Fill(w, fblob, *Q2var);
    Q2_untuned->Fill(*fScaleFactor, fblob, *Q2var);
    if (ist2k) {
      EnuPLepThetaLep->Fill(w, fblob, *Enu_true, *PLep_v, acos(*CosLep));
      EnuPLepThetaLep_untuned->Fill(1, fblob, *Enu_true, *PLep_v,
                                    acos(*CosLep));
    } else {
      EnuPLepEAvHad->Fill(w, fblob, *Enu_true, *PLep_v, *EavAlt);
      EnuQ2EAvHad->Fill(w, fblob, *Enu_true, *Q2var, *EavAlt);
      EnuPtLepEAvHad->Fill(w, fblob, *Enu_true,
                           (*PLep_v) * sqrt(1 - pow(*CosLep, 2)), *EavAlt);
    }
    ent_it++;
  }
}

std::string input_file, output_file;
bool ist2k = false;
bool bymode = false;
int tgta_select = 0;

void SayUsage(char const *argv[]) {
  std::cout << "[USAGE]: " << argv[0]
            << "\n"
               "\t-i <fulldetprocess.root>             : TChain descriptor for"
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
      inpfile = argv[++opt];
    } else if (std::string(argv[opt]) == "-uF") {
      fluxthrowfile = argv[++opt];
    }
  }
}

int main(int argc, char const *argv[]) {
  if (argc < 4) {
    std::cout << "Expects 3 arguments." << std::endl;
    return 1;
  }

  gStyle->SetOptStat(false);

  TFile fin(argv[1]);
  if (fin.IsZombie()) {
    std::cout << "Failed to read " << argv[1] << std::endl;
    return 2;
  }

  TTreeReader rdr("T2KNOvATruthTree", &fin);

  bool ist2k = false;

  if (std::string(argv[2]) == "NOvAND") {
    ist2k = false;
  } else if (std::string(argv[2]) == "ND280") {
    ist2k = true;
  } else {
    std::cout << "[ERROR]: Expected option 2 to specify either NOvAND or ND280."
              << std::endl;
    return 1;
  }

  int tgta_select = 0;
  if (argc > 4) {
    if (std::string(argv[4]) == "C") {
      tgta_select = 12;
    } else if (std::string(argv[4]) == "H") {
      tgta_select = 1;
    } else if (std::string(argv[4]) == "O") {
      tgta_select = 16;
    } else if (std::string(argv[4]) == "any") {
      tgta_select = 0;
    } else {
      std::cout << "Invalid target selector passed: " << argv[4]
                << ". Should be C/H/O/any." << std::endl;
      return 1;
    }
  }

  Fill(rdr, ist2k, tgta_select);

  TFile fout(argv[3], "UPDATE");
  if (fin.IsZombie()) {
    std::cout << "Failed to write " << argv[3] << std::endl;
    return 2;
  }

  TDirectory *dout = &fout;
  if (argc > 5) {
    std::string fqdir = argv[5];
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

  auto axis_labeler = [](TH1D &h) {
    h.GetXaxis()->SetBinLabel(1, "CCInc");
    h.GetXaxis()->SetBinLabel(2, "CC0#pi");
    h.GetXaxis()->SetBinLabel(3, "CC1#pi^{#pm}");
    h.GetXaxis()->SetBinLabel(4, "CC1#pi^{0}");
    h.GetXaxis()->SetBinLabel(5, "CCMulti#pi + CCOther");
    h.GetXaxis()->SetBinLabel(6, "NCInc");
    h.GetXaxis()->SetBinLabel(7, "NC0#pi");
    h.GetXaxis()->SetBinLabel(8, "NC1#pi^{#pm}");
    h.GetXaxis()->SetBinLabel(9, "NC1#pi^{0}");
    h.GetXaxis()->SetBinLabel(10, "NCMulti#pi + CCOther");
  };

  totxsecs->Apply(axis_labeler);
  totxsecs->Write(dout);

  totxsecs_untuned->Apply(axis_labeler);
  totxsecs_untuned->Write(dout);

  Enu->Write(dout);
  Enu_untuned->Write(dout);

  Q2->Write(dout);
  Q2_untuned->Write(dout);

  auto ProjYZ = [=](TH3D const &h) -> TH2D {
    std::unique_ptr<TH2D> h2(dynamic_cast<TH2D *>(h.Project3D("yz")));
    h2->SetDirectory(nullptr);
    h2->SetName(FixProjName(h.GetName(),"yz"));
    return TH2D(*h2.get());
  };
  auto ProjYX = [=](TH3D const &h) -> TH2D {
    std::unique_ptr<TH2D> h2(dynamic_cast<TH2D *>(h.Project3D("yx")));
    h2->SetDirectory(nullptr);
    h2->SetName(FixProjName(h.GetName(),"yx"));
    return TH2D(*h2.get());
  };
  auto ProjXZ = [=](TH3D const &h) -> TH2D {
    std::unique_ptr<TH2D> h2(dynamic_cast<TH2D *>(h.Project3D("xz")));
    h2->SetDirectory(nullptr);
    h2->SetName(FixProjName(h.GetName(),"xz"));
    return TH2D(*h2.get());
  };

  auto ProjX = [=](TH3D const &h) -> TH1D {
    std::unique_ptr<TH1D> h1(dynamic_cast<TH1D *>(h.Project3D("x")));
    h1->SetDirectory(nullptr);
    h1->SetName(FixProjName(h.GetName(),"x"));
    return TH1D(*h1.get());
  };

  auto ProjY = [=](TH3D const &h) -> TH1D {
    std::unique_ptr<TH1D> h1(dynamic_cast<TH1D *>(h.Project3D("y")));
    h1->SetDirectory(nullptr);
    h1->SetName(FixProjName(h.GetName(),"y"));
    return TH1D(*h1.get());
  };

  auto ProjZ = [=](TH3D const &h) -> TH1D {
    std::unique_ptr<TH1D> h1(dynamic_cast<TH1D *>(h.Project3D("z")));
    h1->SetDirectory(nullptr);
    h1->SetName(FixProjName(h.GetName(),"z"));
    return TH1D(*h1.get());
  };

  if (ist2k) {
    EnuPLepThetaLep->Write(dout);
    auto PLepThetaLep = EnuPLepThetaLep->Transform<TH2D>(ProjYZ);
    PLepThetaLep.Write(dout);
    auto EnuPLep = EnuPLepThetaLep->Transform<TH2D>(ProjYX);
    EnuPLep.Write(dout);
    auto EnuThetaLep = EnuPLepThetaLep->Transform<TH2D>(ProjXZ);
    EnuThetaLep.Write(dout);

    auto EnuProj = EnuPLepThetaLep->Transform<TH1D>(ProjX);
    EnuProj.Write(dout);
    auto PLep = EnuPLepThetaLep->Transform<TH1D>(ProjY);
    PLep.Write(dout);
    auto ThetaLep = EnuPLepThetaLep->Transform<TH1D>(ProjZ);
    ThetaLep.Write(dout);

    EnuPLepThetaLep_untuned->Write(dout);
    auto PLepThetaLep_untuned =
        EnuPLepThetaLep_untuned->Transform<TH2D>(ProjYZ);
    PLepThetaLep_untuned.Write(dout);
    auto EnuPLep_untuned = EnuPLepThetaLep_untuned->Transform<TH2D>(ProjYX);
    EnuPLep_untuned.Write(dout);
    auto EnuThetaLep_untuned = EnuPLepThetaLep_untuned->Transform<TH2D>(ProjXZ);
    EnuThetaLep_untuned.Write(dout);

    auto EnuProj_untuned = EnuPLepThetaLep_untuned->Transform<TH1D>(ProjX);
    EnuProj_untuned.Write(dout);
    auto PLep_untuned = EnuPLepThetaLep_untuned->Transform<TH1D>(ProjY);
    PLep_untuned.Write(dout);
    auto ThetaLep_untuned = EnuPLepThetaLep_untuned->Transform<TH1D>(ProjZ);
    ThetaLep_untuned.Write(dout);
  } else {
    EnuPLepEAvHad->Write(dout);
    auto PLepEAvHad = EnuPLepEAvHad->Transform<TH2D>(ProjYZ);
    PLepEAvHad.Write(dout);
    auto EnuPLep = EnuPLepEAvHad->Transform<TH2D>(ProjYX);
    EnuPLep.Write(dout);
    auto EnuEAvHad = EnuPLepEAvHad->Transform<TH2D>(ProjXZ);
    EnuEAvHad.Write(dout);

    EnuQ2EAvHad->Write(dout);
    auto Q2EAvHad = EnuQ2EAvHad->Transform<TH2D>(ProjYZ);
    Q2EAvHad.Write(dout);
    auto EnuQ2 = EnuQ2EAvHad->Transform<TH2D>(ProjYX);
    EnuQ2.Write(dout);

    EnuPtLepEAvHad->Write(dout);
    auto PtLepEAvHad = EnuPtLepEAvHad->Transform<TH2D>(ProjYZ);
    PtLepEAvHad.Write(dout);
    auto EnuPtLep = EnuPtLepEAvHad->Transform<TH2D>(ProjYX);
    EnuPtLep.Write(dout);
  }

  fout.Close();
}
