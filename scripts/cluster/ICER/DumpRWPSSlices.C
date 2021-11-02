#include <iostream>

#include "T2KNOvAROOTHelper.hxx"
#include "TCanvas.h"
#include "TExec.h"
#include "TLatex.h"
#include "TStyle.h"
#include "colordef.h"

int main() {
  gStyle->SetOptStat(false);

  long unsigned int BWRPalette = (long unsigned int)GetBWRPalette();

  std::unique_ptr<TFile> RWHists(new TFile("FDSInputs/FakeDataHists.root"));

  std::unique_ptr<TH3> EnuPLepThetaLep_NEUT_CC0pi =
      GetTH3(RWHists, "NEUT/ND280/C/numu/EnuPLepThetaLep_CC0pi", false);
  std::unique_ptr<TH3> EnuPLepThetaLep_GENIE_CC0pi =
      GetTH3(RWHists, "GENIE/ND280/C/numu/EnuPLepThetaLep_CC0pi", false);

  std::unique_ptr<TH3> EnuPLepThetaLep_RW_CC0pi(static_cast<TH3*>(
      EnuPLepThetaLep_NEUT_CC0pi->Clone("EnuPLepThetaLep_RW_CC0pi")));
  EnuPLepThetaLep_RW_CC0pi->Divide(EnuPLepThetaLep_GENIE_CC0pi.get());

  std::unique_ptr<TH3> EnuPLepThetaLep_RW_CC0pi_scrub(static_cast<TH3*>(
      EnuPLepThetaLep_RW_CC0pi->Clone("EnuPLepThetaLep_RW_CC0pi_scrub")));
  ScrubLowStatsBins(EnuPLepThetaLep_NEUT_CC0pi, EnuPLepThetaLep_GENIE_CC0pi,
                    EnuPLepThetaLep_RW_CC0pi_scrub, 1.0 / sqrt(25));

  TCanvas c1("c1", "", 1000, 1000);
  c1.Print("CC0PiPS.pdf[");

  for (int i = 0; i < EnuPLepThetaLep_NEUT_CC0pi->GetXaxis()->GetNbins(); ++i) {
    c1.Clear();

    TPad p1("p1", "", 0, 0.6, 0.5, 0.9);
    p1.AppendPad();

    c1.cd();
    TPad p2("p2", "", 0.5, 0.6, 1, 0.9);
    p2.AppendPad();

    c1.cd();
    TPad p3("p3", "", 0, 0.3, 0.5, 0.6);
    p3.AppendPad();

    c1.cd();
    TPad p4("p4", "", 0.5, 0.3, 1, 0.6);
    p4.AppendPad();

    c1.cd();
    TPad p5("p5", "", 0, 0, 0.5, 0.3);
    p5.AppendPad();

    c1.cd();
    TPad p6("p6", "", 0.5, 0, 1, 0.3);
    p6.AppendPad();

    std::unique_ptr<TH2> PLepThetaLep_NEUT_CC0pi =
        ProjectYZSlice(EnuPLepThetaLep_NEUT_CC0pi, i + 1, i + 2);
    PLepThetaLep_NEUT_CC0pi->SetName("PLepThetaLep_NEUT_CC0pi");
    std::unique_ptr<TH2> PLepThetaLep_GENIE_CC0pi =
        ProjectYZSlice(EnuPLepThetaLep_GENIE_CC0pi, i + 1, i + 2);
    PLepThetaLep_GENIE_CC0pi->SetName("PLepThetaLep_GENIE_CC0pi");

    std::unique_ptr<TH2> PLepThetaLep_RW_CC0pi =
        ProjectYZSlice(EnuPLepThetaLep_RW_CC0pi, i + 1, i + 2);
    PLepThetaLep_RW_CC0pi->SetName("PLepThetaLep_RW_CC0pi");

    std::unique_ptr<TH2> PLepThetaLep_RW_CC0pi_scrub =
        ProjectYZSlice(EnuPLepThetaLep_RW_CC0pi_scrub, i + 1, i + 2);
    PLepThetaLep_RW_CC0pi_scrub->SetName("PLepThetaLep_RW_CC0pi_scrub");

    double max = std::max(PLepThetaLep_NEUT_CC0pi->GetMaximum(),
                          PLepThetaLep_GENIE_CC0pi->GetMaximum());

    PLepThetaLep_NEUT_CC0pi->GetZaxis()->SetRangeUser(0, max);
    PLepThetaLep_GENIE_CC0pi->GetZaxis()->SetRangeUser(0, max);

    std::unique_ptr<TH2> NEUT_on_GENIE =
        GetNonOverlap2D(PLepThetaLep_NEUT_CC0pi, PLepThetaLep_GENIE_CC0pi);
    NEUT_on_GENIE->SetName("NEUT_on_GENIE");
    std::unique_ptr<TH2> GENIE_ON_NEUT =
        GetNonOverlap2D(PLepThetaLep_GENIE_CC0pi, PLepThetaLep_NEUT_CC0pi);
    GENIE_ON_NEUT->SetName("GENIE_ON_NEUT");

    NEUT_on_GENIE->GetZaxis()->SetRangeUser(0, max);
    GENIE_ON_NEUT->GetZaxis()->SetRangeUser(0, max);

    std::string BWR_cmd = "gStyle->SetPalette(100, (int*)" +
                          std::to_string(BWRPalette) +
                          ");gStyle->SetNumberContours(100);";

    TExec* ex1 = new TExec("ex1", "gStyle->SetPalette(kBird);");

    p1.cd();
    PLepThetaLep_NEUT_CC0pi->SetTitle("");
    auto cPLepThetaLep_NEUT_CC0pi = PLepThetaLep_NEUT_CC0pi->DrawClone("COLZ");
    ex1->Draw();
    cPLepThetaLep_NEUT_CC0pi->DrawClone("COLZ SAME");

    p2.cd();
    PLepThetaLep_GENIE_CC0pi->SetTitle("");
    auto cPLepThetaLep_GENIE_CC0pi =
        PLepThetaLep_GENIE_CC0pi->DrawClone("COLZ");
    ex1->Draw();
    cPLepThetaLep_GENIE_CC0pi->DrawClone("COLZ SAME");

    p3.cd();
    PLepThetaLep_NEUT_CC0pi->SetLineColor(kRed);
    PLepThetaLep_GENIE_CC0pi->SetLineColor(kBlue);
    PLepThetaLep_GENIE_CC0pi->SetLineStyle(3);
    PLepThetaLep_NEUT_CC0pi->DrawClone("BOX");
    PLepThetaLep_GENIE_CC0pi->DrawClone("BOXSAME");

    p4.cd();
    NEUT_on_GENIE->SetLineColor(kRed);
    GENIE_ON_NEUT->SetLineColor(kBlue);
    GENIE_ON_NEUT->SetLineStyle(3);
    NEUT_on_GENIE->DrawClone("BOX");
    GENIE_ON_NEUT->DrawClone("BOXSAME");

    p5.cd();
    TExec* ex2 = new TExec("ex2", BWR_cmd.c_str());
    PLepThetaLep_RW_CC0pi->GetZaxis()->SetRangeUser(-5, 5);
    PLepThetaLep_RW_CC0pi->SetTitle("");
    auto cPLepThetaLep_RW_CC0pi = PLepThetaLep_RW_CC0pi->DrawClone("COLZ");
    ex2->Draw();
    cPLepThetaLep_RW_CC0pi->DrawClone("COLZ SAME");

    p6.cd();
    PLepThetaLep_RW_CC0pi_scrub->GetZaxis()->SetRangeUser(-5, 5);
    PLepThetaLep_RW_CC0pi_scrub->SetTitle("");
    auto cPLepThetaLep_RW_CC0pi_scrub =
        PLepThetaLep_RW_CC0pi_scrub->DrawClone("COLZ");
    ex2->Draw();
    cPLepThetaLep_RW_CC0pi_scrub->DrawClone("COLZ SAME");

    c1.cd();
    TLatex ltx;
    ltx.SetTextAlign(22);
    ltx.DrawLatexNDC(
        0.5, 0.95,
        (std::to_string(
             EnuPLepThetaLep_NEUT_CC0pi->GetXaxis()->GetBinLowEdge(i + 1)) +
         " < Enu < " +
         std::to_string(
             EnuPLepThetaLep_NEUT_CC0pi->GetXaxis()->GetBinUpEdge(i + 1)) +
         " GeV")
            .c_str());

    ltx.SetTextColor(kRed);
    ltx.DrawLatexNDC(0.25, 0.8, "NEUT");

    ltx.SetTextColor(kBlue);
    ltx.DrawLatexNDC(0.75, 0.8, "GENIE");

    ltx.SetTextColor(kBlack);
    ltx.SetTextSize(0.03);
    ltx.DrawLatexNDC(0.75, 0.5, "Exclusive phase space");

    ltx.DrawLatexNDC(0.25, 0.2, "Weighting");
    ltx.DrawLatexNDC(0.75, 0.2, "Weighting (25 ev/bin cut)");

    c1.Print("CC0PiPS.pdf");
  }

  c1.Print("CC0PiPS.pdf]");
}