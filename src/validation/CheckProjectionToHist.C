#include "T2KNOvA/FakeDataHelper.hxx"
#include "T2KNOvA/ROOTHelper.hxx"
#include "T2KNOvA/TrueSelectionHelper.hxx"

#include "colordef.h"
#include "plotutils.h"

#include "TCanvas.h"
#include "TFile.h"
#include "TH3.h"
#include "TStyle.h"

#include <memory>

int main(int argc, char const *argv[]) {
  gStyle->SetOptStat(false);

  long unsigned int BWRPalette = (long unsigned int)GetBWRPalette();

  std::unique_ptr<TFile> RWHists(new TFile(argv[1]));
  std::unique_ptr<TFile> ValidHists(new TFile(argv[2]));

  TCanvas c1("c1", "", 1200, 800);
  c1.Print("PLepProjections.pdf[");

  for (int sel : {
           t2knova::kCCInc,
           t2knova::kCC0pi,
           t2knova::kCC1cpi,
           t2knova::kCC1pi0,
           t2knova::kCCmultipi,
           t2knova::kCC1Gamma,
           t2knova::kCCOther,
       }) {
    std::unique_ptr<TH3> EnuPLepThetaLep_NEUT = GetTH<TH3>(
        RWHists,
        "NEUT/ND280/C/numu/EnuPLepThetaLep_" + t2knova::SelectionList[sel],
        false);
    std::unique_ptr<TH3> EnuPLepThetaLep_GENIE = GetTH<TH3>(
        RWHists,
        "GENIE/ND280/C/numu/EnuPLepThetaLep_" + t2knova::SelectionList[sel],
        false);

    std::unique_ptr<TH1> PLep_NEUT = GetTH<TH1>(
        RWHists, "NEUT/ND280/C/numu/PLep_" + t2knova::SelectionList[sel],
        false);
    std::unique_ptr<TH1> PLep_GENIE = GetTH<TH1>(
        RWHists, "GENIE/ND280/C/numu/PLep_" + t2knova::SelectionList[sel],
        false);

    std::unique_ptr<TH3> EnuPLepThetaLep_NEUT_to_GENIE_t2k_Valid =
        GetTH<TH3>(ValidHists,
                   "ND280/T2KND_to_NOvA/C/numu/EnuPLepThetaLep_" +
                       t2knova::SelectionList[sel],
                   false);

    std::unique_ptr<TH3> EnuPLepThetaLep_NEUT_t2k_Valid = GetTH<TH3>(
        ValidHists,
        "ND280/T2KNDTune/C/numu/EnuPLepThetaLep_" + t2knova::SelectionList[sel],
        false);
    std::unique_ptr<TH3> EnuPLepThetaLep_GENIE_t2k_Valid = GetTH<TH3>(
        ValidHists,
        "ND280/NOvATune/C/numu/EnuPLepThetaLep_" + t2knova::SelectionList[sel],
        false);

    std::unique_ptr<TH1> PLep_NEUT_t2k_Valid = GetTH<TH1>(
        ValidHists,
        "ND280/T2KNDTune/C/numu/PLep_" + t2knova::SelectionList[sel], false);
    std::unique_ptr<TH1> PLep_GENIE_t2k_Valid = GetTH<TH1>(
        ValidHists, "ND280/NOvATune/C/numu/PLep_" + t2knova::SelectionList[sel],
        false);

    c1.cd();

    std::unique_ptr<TH1> PLep_proj_NEUT = ProjectY_XSlice(
        EnuPLepThetaLep_NEUT, 1, EnuPLepThetaLep_NEUT->GetXaxis()->GetNbins());
    PLep_proj_NEUT->SetName("PLep_proj_NEUT");

    std::unique_ptr<TH1> PLep_proj_GENIE = ProjectY_XSlice(
        EnuPLepThetaLep_GENIE, 1, EnuPLepThetaLep_NEUT->GetXaxis()->GetNbins());
    PLep_proj_GENIE->SetName("PLep_proj_GENIE");

    PLep_proj_NEUT->SetTitle(
        ("PLep Projection, T2K Flux, " + t2knova::SelectionList[sel]).c_str());

    StyleAxes(PLep_proj_NEUT);
    StyleTH1Line(PLep_proj_NEUT, kBlue, 2, 1);
    StyleTH1Line(PLep_proj_GENIE, kRed, 2, 1);
    StyleTH1Line(PLep_NEUT, kBlue, 2, 2);
    StyleTH1Line(PLep_GENIE, kRed, 2, 2);

    DrawTH1s({PLep_proj_NEUT, PLep_proj_GENIE, PLep_NEUT, PLep_GENIE}, "HIST");

    // std::unique_ptr<TH1> PLep_proj_NEUT_t2k_Valid =
    //     ProjectY_XSlice(EnuPLepThetaLep_NEUT_t2k_Valid, 1,
    //                     EnuPLepThetaLep_NEUT->GetXaxis()->GetNbins());
    // PLep_proj_NEUT_t2k_Valid->SetName("PLep_proj_NEUT_t2k_Valid");

    // std::unique_ptr<TH1> PLep_proj_GENIE_t2k_Valid =
    //     ProjectY_XSlice(EnuPLepThetaLep_GENIE_t2k_Valid, 1,
    //                     EnuPLepThetaLep_NEUT->GetXaxis()->GetNbins());
    // PLep_proj_GENIE_t2k_Valid->SetName("PLep_proj_GENIE_t2k_Valid");

    // PLep_proj_NEUT_t2k_Valid->SetTitle(
    //     ("PLep Projection, T2K Flux, " +
    //     t2knova::SelectionList[sel]).c_str());

    // StyleAxes(PLep_proj_NEUT_t2k_Valid);
    // StyleTH1Line(PLep_proj_NEUT_t2k_Valid, kBlue, 2, 1);
    // StyleTH1Line(PLep_proj_GENIE_t2k_Valid, kRed, 2, 1);
    // StyleTH1Line(PLep_NEUT_t2k_Valid, kBlue, 2, 2);
    // StyleTH1Line(PLep_GENIE_t2k_Valid, kRed, 2, 2);

    // DrawTH1s({PLep_proj_NEUT_t2k_Valid, PLep_proj_GENIE_t2k_Valid,
    //           PLep_NEUT_t2k_Valid, PLep_GENIE_t2k_Valid},
    //          "HIST");

    c1.Print("PLepProjections.pdf");
  }
  c1.Print("PLepProjections.pdf]");
}