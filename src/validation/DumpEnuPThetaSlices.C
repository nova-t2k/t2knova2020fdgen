#include <iostream>

#include "T2KNOvA/FakeDataHelper.hxx"
#include "T2KNOvA/ROOTHelper.hxx"
#include "T2KNOvA/TrueSelectionHelper.hxx"

#include "colordef.h"
#include "plotutils.h"

#include "TCanvas.h"
#include "TExec.h"
#include "TLatex.h"
#include "TStyle.h"

double low_z = 0;
double up_z = 2;

double upplep = 3;

int main(int argc, char const *argv[]) {
  gStyle->SetOptStat(false);

  long unsigned int BWRPalette = (long unsigned int)GetBWRPalette();

  std::unique_ptr<TFile> RWHists(new TFile(argv[1]));
  std::unique_ptr<TFile> ValidHists(new TFile(argv[2]));

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

    std::unique_ptr<TH3> EnuPLepThetaLep_RW(
        static_cast<TH3 *>(EnuPLepThetaLep_GENIE->Clone("EnuPLepThetaLep_RW")));
    EnuPLepThetaLep_RW->Divide(EnuPLepThetaLep_NEUT.get());

    std::unique_ptr<TH3> EnuPLepThetaLep_RW_scrub(static_cast<TH3 *>(
        EnuPLepThetaLep_RW->Clone("EnuPLepThetaLep_RW_scrub")));
    ScrubLowStatsBins(EnuPLepThetaLep_NEUT, EnuPLepThetaLep_GENIE,
                      EnuPLepThetaLep_RW_scrub, t2knova::MaxFracError);

    TCanvas c1("c1", "", 1000, 1000);
    c1.Print(("PS_" + t2knova::SelectionList[sel] + ".pdf[").c_str());

    for (int i = 0; i < EnuPLepThetaLep_NEUT->GetXaxis()->GetNbins(); ++i) {
      c1.Clear();

      std::unique_ptr<TH2> PLepThetaLep_NEUT =
          ProjectYZSlice(EnuPLepThetaLep_NEUT, i + 1, i + 2);
      PLepThetaLep_NEUT->SetName("PLepThetaLep_NEUT");
      std::unique_ptr<TH2> PLepThetaLep_GENIE =
          ProjectYZSlice(EnuPLepThetaLep_GENIE, i + 1, i + 2);
      PLepThetaLep_GENIE->SetName("PLepThetaLep_GENIE");

      std::unique_ptr<TH2> PLepThetaLep_RW =
          ProjectYZSlice(EnuPLepThetaLep_RW, i + 1, i + 2);
      PLepThetaLep_RW->SetName("PLepThetaLep_RW");

      std::unique_ptr<TH2> PLepThetaLep_RW_scrub =
          ProjectYZSlice(EnuPLepThetaLep_RW_scrub, i + 1, i + 2);
      PLepThetaLep_RW_scrub->SetName("PLepThetaLep_RW_scrub");

      double max = std::max(PLepThetaLep_NEUT->GetMaximum(),
                            PLepThetaLep_GENIE->GetMaximum());

      std::string BWR_cmd = "gStyle->SetPalette(100, (int*)" +
                            std::to_string(BWRPalette) +
                            ");gStyle->SetNumberContours(100);";

      TExec *ex1 = new TExec("ex1", "gStyle->SetPalette(kBird);");
      TExec *ex2 = new TExec("ex2", BWR_cmd.c_str());

      {
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

        PLepThetaLep_NEUT->GetZaxis()->SetRangeUser(0, max);
        PLepThetaLep_GENIE->GetZaxis()->SetRangeUser(0, max);

        std::unique_ptr<TH2> NEUT_on_GENIE =
            GetNonOverlap2D(PLepThetaLep_NEUT, PLepThetaLep_GENIE);
        NEUT_on_GENIE->SetName("NEUT_on_GENIE");
        std::unique_ptr<TH2> GENIE_ON_NEUT =
            GetNonOverlap2D(PLepThetaLep_GENIE, PLepThetaLep_NEUT);
        GENIE_ON_NEUT->SetName("GENIE_ON_NEUT");

        NEUT_on_GENIE->GetZaxis()->SetRangeUser(0, max);
        GENIE_ON_NEUT->GetZaxis()->SetRangeUser(0, max);

        PLepThetaLep_NEUT->GetYaxis()->SetRangeUser(0, upplep);
        PLepThetaLep_GENIE->GetYaxis()->SetRangeUser(0, upplep);
        NEUT_on_GENIE->GetYaxis()->SetRangeUser(0, upplep);
        PLepThetaLep_RW->GetYaxis()->SetRangeUser(0, upplep);
        PLepThetaLep_RW_scrub->GetYaxis()->SetRangeUser(0, upplep);

        p1.cd();
        PLepThetaLep_NEUT->SetTitle("");
        auto cPLepThetaLep_NEUT = PLepThetaLep_NEUT->DrawClone("COLZ");
        ex1->Draw();
        cPLepThetaLep_NEUT->DrawClone("COLZ SAME");

        p2.cd();
        PLepThetaLep_GENIE->SetTitle("");
        auto cPLepThetaLep_GENIE = PLepThetaLep_GENIE->DrawClone("COLZ");
        ex1->Draw();
        cPLepThetaLep_GENIE->DrawClone("COLZ SAME");

        p3.cd();
        PLepThetaLep_NEUT->SetLineColor(kRed);
        PLepThetaLep_GENIE->SetLineColor(kBlue);
        PLepThetaLep_GENIE->SetLineStyle(3);
        PLepThetaLep_NEUT->DrawClone("BOX");
        PLepThetaLep_GENIE->DrawClone("BOXSAME");

        p4.cd();
        NEUT_on_GENIE->SetTitle("");
        NEUT_on_GENIE->SetLineColor(kRed);
        GENIE_ON_NEUT->SetLineColor(kBlue);
        GENIE_ON_NEUT->SetLineStyle(3);
        NEUT_on_GENIE->DrawClone("BOX");
        GENIE_ON_NEUT->DrawClone("BOXSAME");

        p5.cd();

        PLepThetaLep_RW->GetZaxis()->SetRangeUser(low_z, up_z);
        PLepThetaLep_RW->SetTitle("");
        auto cPLepThetaLep_RW = PLepThetaLep_RW->DrawClone("COLZ");
        ex2->Draw();
        cPLepThetaLep_RW->DrawClone("COLZ SAME");

        p6.cd();
        PLepThetaLep_RW_scrub->GetZaxis()->SetRangeUser(low_z, up_z);
        PLepThetaLep_RW_scrub->SetTitle("");
        auto cPLepThetaLep_RW_scrub = PLepThetaLep_RW_scrub->DrawClone("COLZ");
        ex2->Draw();
        cPLepThetaLep_RW_scrub->DrawClone("COLZ SAME");

        c1.cd();
        TLatex ltx;
        ltx.SetTextAlign(22);
        ltx.DrawLatexNDC(
            0.5, 0.95,
            (std::to_string(
                 EnuPLepThetaLep_NEUT->GetXaxis()->GetBinLowEdge(i + 1)) +
             " < Enu < " +
             std::to_string(
                 EnuPLepThetaLep_NEUT->GetXaxis()->GetBinUpEdge(i + 1)) +
             " GeV")
                .c_str());

        ltx.SetTextColor(kRed);
        ltx.DrawLatexNDC(0.25, 0.9, "NEUT");

        ltx.SetTextColor(kBlue);
        ltx.DrawLatexNDC(0.75, 0.9, "GENIE");

        ltx.SetTextColor(kBlack);
        ltx.SetTextSize(0.03);
        ltx.DrawLatexNDC(0.75, 0.59, "Exclusive phase space");

        ltx.SetTextAlign(12);
        ltx.SetTextSize(0.015);
        ltx.DrawLatexNDC(0.1, 0.6,
                         (std::string("NEUT Integral: ") +
                          std::to_string(IntegralTH2(PLepThetaLep_NEUT,
                                                     t2knova::MaxFracError)) +
                          "/" + std::to_string(IntegralTH2(PLepThetaLep_NEUT)))
                             .c_str());
        ltx.DrawLatexNDC(0.1, 0.58,
                         (std::string("GENIE Integral: ") +
                          std::to_string(IntegralTH2(PLepThetaLep_GENIE,
                                                     t2knova::MaxFracError)) +
                          "/" + std::to_string(IntegralTH2(PLepThetaLep_GENIE)))
                             .c_str());

        ltx.SetTextAlign(22);
        ltx.SetTextSize(0.03);
        ltx.DrawLatexNDC(0.25, 0.3, "Weighting");
        ltx.DrawLatexNDC(0.75, 0.3,
                         (std::string("Weighting (") +
                          std::to_string(int(t2knova::NMinEvs)) +
                          " ev/bin cut)")
                             .c_str());

        c1.Print(("PS_" + t2knova::SelectionList[sel] + ".pdf").c_str());
      }

      //////////////////////////////////////

      {
        gStyle->SetPalette(kBird);

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
        std::unique_ptr<TH2> PLepThetaLep_NEUT_to_GENIE_t2k_Valid =
            ProjectYZSlice(EnuPLepThetaLep_NEUT_to_GENIE_t2k_Valid, i + 1,
                           i + 2);
        PLepThetaLep_NEUT_to_GENIE_t2k_Valid->SetName(
            "PLepThetaLep_NEUT_to_GENIE_t2k_Valid");

        std::unique_ptr<TH2> PLepThetaLep_NEUT_t2k_Valid =
            ProjectYZSlice(EnuPLepThetaLep_NEUT_t2k_Valid, i + 1, i + 2);
        PLepThetaLep_NEUT_t2k_Valid->SetName("PLepThetaLep_NEUT_t2k_Valid");

        std::unique_ptr<TH2> PLepThetaLep_GENIE_t2k_Valid =
            ProjectYZSlice(EnuPLepThetaLep_GENIE_t2k_Valid, i + 1, i + 2);
        PLepThetaLep_GENIE_t2k_Valid->SetName("PLepThetaLep_GENIE_t2k_Valid");

        std::unique_ptr<TH1> PLep_NEUT_to_GENIE_t2k_Valid = ProjectY_XSlice(
            EnuPLepThetaLep_NEUT_to_GENIE_t2k_Valid, i + 1, i + 2);
        PLep_NEUT_to_GENIE_t2k_Valid->SetName("PLep_NEUT_to_GENIE_t2k_Valid");

        std::unique_ptr<TH1> PLep_NEUT_t2k_Valid =
            ProjectY_XSlice(EnuPLepThetaLep_NEUT_t2k_Valid, i + 1, i + 2);
        PLep_NEUT_t2k_Valid->SetName("PLep_NEUT_t2k_Valid");

        std::unique_ptr<TH1> PLep_GENIE_t2k_Valid =
            ProjectY_XSlice(EnuPLepThetaLep_GENIE_t2k_Valid, i + 1, i + 2);
        PLep_GENIE_t2k_Valid->SetName("PLep_GENIE_t2k_Valid");

        std::unique_ptr<TH1> ThetaLep_NEUT_to_GENIE_t2k_Valid = ProjectZ_XSlice(
            EnuPLepThetaLep_NEUT_to_GENIE_t2k_Valid, i + 1, i + 2);
        ThetaLep_NEUT_to_GENIE_t2k_Valid->SetName(
            "ThetaLep_NEUT_to_GENIE_t2k_Valid");

        std::unique_ptr<TH1> ThetaLep_NEUT_t2k_Valid =
            ProjectZ_XSlice(EnuPLepThetaLep_NEUT_t2k_Valid, i + 1, i + 2);
        ThetaLep_NEUT_t2k_Valid->SetName("ThetaLep_NEUT_t2k_Valid");

        std::unique_ptr<TH1> ThetaLep_GENIE_t2k_Valid =
            ProjectZ_XSlice(EnuPLepThetaLep_GENIE_t2k_Valid, i + 1, i + 2);
        ThetaLep_GENIE_t2k_Valid->SetName("ThetaLep_GENIE_t2k_Valid");

        PLepThetaLep_NEUT_t2k_Valid->GetZaxis()->SetRangeUser(0, max);
        PLepThetaLep_GENIE_t2k_Valid->GetZaxis()->SetRangeUser(0, max);

        PLepThetaLep_NEUT->GetYaxis()->SetRangeUser(0, upplep);
        PLepThetaLep_GENIE->GetYaxis()->SetRangeUser(0, upplep);
        PLepThetaLep_GENIE_t2k_Valid->GetYaxis()->SetRangeUser(0, upplep);
        PLepThetaLep_RW_scrub->GetYaxis()->SetRangeUser(0, upplep);
        PLepThetaLep_NEUT_to_GENIE_t2k_Valid->GetYaxis()->SetRangeUser(0,
                                                                       upplep);

        double int_NEUT_TO_GENIE =
            IntegralTH(PLepThetaLep_NEUT_to_GENIE_t2k_Valid);
        double int_GENIE = IntegralTH(PLepThetaLep_GENIE_t2k_Valid);
        double int_NEUT = IntegralTH(PLepThetaLep_NEUT_t2k_Valid);

        p1.cd();
        ThetaLep_NEUT_t2k_Valid->SetLineColor(kRed);
        ThetaLep_GENIE_t2k_Valid->SetLineColor(kBlue);
        ThetaLep_NEUT_t2k_Valid->SetLineWidth(2);
        ThetaLep_GENIE_t2k_Valid->SetLineWidth(2);
        ThetaLep_NEUT_t2k_Valid->SetTitle("");
        ThetaLep_NEUT_to_GENIE_t2k_Valid->SetLineColor(kGreen);
        DrawTH1s(
            std::vector<std::reference_wrapper<std::unique_ptr<TH1>>>{
                ThetaLep_NEUT_t2k_Valid, ThetaLep_GENIE_t2k_Valid,
                ThetaLep_NEUT_to_GENIE_t2k_Valid},
            "EHIST");

        p3.cd();
        std::unique_ptr<TH2> PLepThetaLep_NEUT_to_GENIE_t2k_Valid_c(
            static_cast<TH2 *>(PLepThetaLep_NEUT_to_GENIE_t2k_Valid->Clone(
                "PLepThetaLep_NEUT_to_GENIE_t2k_Valid_c")));
        PLepThetaLep_NEUT_to_GENIE_t2k_Valid_c->Divide(
            PLepThetaLep_GENIE.get());
        PLepThetaLep_NEUT_to_GENIE_t2k_Valid_c->GetZaxis()->SetRangeUser(0.8,
                                                                         1.2);
        PLepThetaLep_NEUT_to_GENIE_t2k_Valid_c->SetTitle("");
        PLepThetaLep_NEUT_to_GENIE_t2k_Valid_c->DrawClone("COLZ");
        ex2->Draw();
        PLepThetaLep_NEUT_to_GENIE_t2k_Valid_c->DrawClone("COLZ SAME");

        p2.cd();
        PLep_NEUT_t2k_Valid->SetLineColor(kRed);
        PLep_GENIE_t2k_Valid->SetLineColor(kBlue);
        PLep_NEUT_t2k_Valid->SetLineWidth(2);
        PLep_GENIE_t2k_Valid->SetLineWidth(2);
        PLep_NEUT_t2k_Valid->SetTitle("");
        PLep_NEUT_to_GENIE_t2k_Valid->SetLineColor(kGreen);
        DrawTH1s(
            std::vector<std::reference_wrapper<std::unique_ptr<TH1>>>{
                PLep_NEUT_t2k_Valid, PLep_GENIE_t2k_Valid,
                PLep_NEUT_to_GENIE_t2k_Valid},
            "EHIST");

        p4.cd();
        PLepThetaLep_GENIE_t2k_Valid->Divide(PLepThetaLep_NEUT_t2k_Valid.get());
        PLepThetaLep_GENIE_t2k_Valid->GetZaxis()->SetRangeUser(low_z, up_z);
        PLepThetaLep_GENIE_t2k_Valid->SetTitle("");
        PLepThetaLep_GENIE_t2k_Valid->DrawClone("COLZ");
        ex2->Draw();
        PLepThetaLep_GENIE_t2k_Valid->DrawClone("COLZ SAME");

        p5.cd();
        PLepThetaLep_RW_scrub->DrawClone("COLZ");
        ex2->Draw();
        PLepThetaLep_RW_scrub->DrawClone("COLZ SAME");

        p6.cd();
        PLepThetaLep_NEUT_to_GENIE_t2k_Valid->Divide(PLepThetaLep_NEUT.get());
        PLepThetaLep_NEUT_to_GENIE_t2k_Valid->GetZaxis()->SetRangeUser(low_z,
                                                                       up_z);
        PLepThetaLep_NEUT_to_GENIE_t2k_Valid->SetTitle("");
        PLepThetaLep_NEUT_to_GENIE_t2k_Valid->DrawClone("COLZ");
        ex2->Draw();
        PLepThetaLep_NEUT_to_GENIE_t2k_Valid->DrawClone("COLZ SAME");

        c1.cd();
        TLatex ltx;
        ltx.SetTextAlign(22);
        ltx.DrawLatexNDC(
            0.5, 0.95,
            (std::to_string(
                 EnuPLepThetaLep_NEUT->GetXaxis()->GetBinLowEdge(i + 1)) +
             " < Enu < " +
             std::to_string(
                 EnuPLepThetaLep_NEUT->GetXaxis()->GetBinUpEdge(i + 1)) +
             " GeV")
                .c_str());

        ltx.SetTextSize(0.03);

        ltx.DrawLatexNDC(0.75, 0.9, "PLep Projection");
        ltx.DrawLatexNDC(0.25, 0.9, "ThetaLep Projection");

        ltx.SetTextAlign(31);
        ltx.SetTextSize(0.02);
        ltx.SetTextColor(kRed);
        ltx.DrawLatexNDC(0.94, 0.75, "NEUT");
        ltx.SetTextColor(kBlue);
        ltx.DrawLatexNDC(0.94, 0.725, "GENIE");
        ltx.SetTextColor(kGreen);
        ltx.DrawLatexNDC(0.94, 0.7, "NEUT->GENIE");

        ltx.SetTextColor(kBlack);
        ltx.SetTextAlign(22);
        ltx.DrawLatexNDC(0.25, 0.6, "EvWeighted/GENIE");
        ltx.DrawLatexNDC(0.75, 0.6, "GENIE/NEUT (Valid)");

        ltx.DrawLatexNDC(0.25, 0.3, "GENIE/NEUT (Inputs, Scrubbed)");
        ltx.DrawLatexNDC(0.75, 0.3, "EvWeighted/NEUT");

        ltx.SetTextAlign(31);

        ltx.DrawLatexNDC(
            0.94, 0.54,
            (std::string("NEUT Integral: ") + std::to_string(int_NEUT))
                .c_str());
        ltx.DrawLatexNDC(
            0.94, 0.52,
            (std::string("GENIE Integral: ") + std::to_string(int_GENIE))
                .c_str());
        ltx.DrawLatexNDC(
            0.94, 0.5,
            (std::string("N->G Integral: ") + std::to_string(int_NEUT_TO_GENIE))
                .c_str());

        c1.Print(("PS_" + t2knova::SelectionList[sel] + ".pdf").c_str());
      }

      std::cout << "Bin " << i << "/"
                << EnuPLepThetaLep_NEUT->GetXaxis()->GetNbins() << std::endl;
    }

    c1.Print(("PS_" + t2knova::SelectionList[sel] + ".pdf]").c_str());
  }
}