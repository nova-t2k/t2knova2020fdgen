#include "T2KNOvAFakeDataHelper.hxx"

#include "colordef.h"
#include "plotutils.h"

#include "TLatex.h"

int main() {
  DeclareColors();

  TFile F_FDH("FakeDataHists.root", "READ");
  TFile F_FDI("FakeDataInputs.root", "READ");

  TCanvas *c1 = MakeCanvasTopLegend();
  c1->Print("T2KNOvATunePreds.pdf[");

  for (auto detstr : {"ND280"}) {
    for (auto selstr :
         {"CCInc", "CC0Pi", "CC1CPi", "CC1Pi0", "CCOther", "totxsecs"}) {
      for (auto projstr : {"Enu", "Q2"}) {
        for (auto tgtstr : {"C", "O", "H"}) {
          for (auto nuspec : {t2knova::kNuMu, t2knova::kNuMub, t2knova::kNuE,
                              t2knova::kNuEb}) {

            if ((std::string(tgtstr) == "H") &&
                ((nuspec == t2knova::kNuMu) || (nuspec == t2knova::kNuE))) {
              continue;
            }

            if (selstr == "totxsecs") {
              TH1 *NEUT = GetTH1(
                  &F_FDH, std::string("NEUT/") + detstr + "/" + tgtstr + "/" +
                              t2knova::all_nuspecies[nuspec] + "/" + selstr);
              if (!NEUT) {
                continue;
              }

              TH1 *NEUT_untuned =
                  GetTH1(&F_FDH, std::string("NEUT/") + detstr + "/" + tgtstr +
                                     "/" + t2knova::all_nuspecies[nuspec] +
                                     "/" + selstr + "_untuned");

              TH1 *GENIE = GetTH1(
                  &F_FDH, std::string("GENIE/") + detstr + "/" + tgtstr + "/" +
                              t2knova::all_nuspecies[nuspec] + "/" + selstr);

              TH1 *GENIE_untuned =
                  GetTH1(&F_FDH, std::string("GENIE/") + detstr + "/" + tgtstr +
                                     "/" + t2knova::all_nuspecies[nuspec] +
                                     "/" + selstr + "_untuned");

              c1->Clear();
              NEUT->GetYaxis()->SetTitle("#sigma_{tot.}");
              StyleTH1Line(NEUT, kT2KRed, 2, 1);
              StyleTH1Line(NEUT_untuned, kT2KRed, 2, 2);
              StyleTH1Line(GENIE, kNOvABlue, 2, 1);
              StyleTH1Line(GENIE_untuned, kNOvABlue, 2, 2);

              auto *leg = MakeTopLegend();
              leg->SetTextSize(0.04);

              leg->AddEntry(NEUT, "NEUT+BANFF", "l");
              leg->AddEntry(NEUT_untuned, "NEUT", "l");
              leg->AddEntry(GENIE, "GENIE+NOvA2020", "l");
              leg->AddEntry(GENIE_untuned, "GENIE", "l");

              DrawTH1s({NEUT, NEUT_untuned, GENIE, GENIE_untuned}, "EHIST");
              leg->Draw();

              TLatex ltx;
              ltx.SetTextAlign(31);
              ltx.DrawLatexNDC(0.9, 0.75,
                               (std::string(detstr) + " " + selstr + " " +
                                tgtstr + " " + t2knova::all_nuspecies[nuspec])
                                   .c_str());

              c1->Print("T2KNOvATunePreds.pdf");
            } else {

              TH1 *NEUT =
                  GetTH1(&F_FDH, std::string("NEUT/") + detstr + "/" + tgtstr +
                                     "/" + t2knova::all_nuspecies[nuspec] +
                                     "/" + projstr + "_" + selstr);
              if (!NEUT) {
                continue;
              }

              TH1 *NEUT_untuned =
                  GetTH1(&F_FDH, std::string("NEUT/") + detstr + "/" + tgtstr +
                                     "/" + t2knova::all_nuspecies[nuspec] +
                                     "/" + projstr + "_untuned_" + selstr);

              TH1 *GENIE =
                  GetTH1(&F_FDH, std::string("GENIE/") + detstr + "/" + tgtstr +
                                     "/" + t2knova::all_nuspecies[nuspec] +
                                     "/" + projstr + "_" + selstr);

              TH1 *GENIE_untuned =
                  GetTH1(&F_FDH, std::string("GENIE/") + detstr + "/" + tgtstr +
                                     "/" + t2knova::all_nuspecies[nuspec] +
                                     "/" + projstr + "_untuned_" + selstr);

              TH1 *Inputs1D = GetTH1(
                  &F_FDI, std::string("t2knd_to_nova_") + projstr + "_" +
                              tgtstr + "_" + t2knova::all_nuspecies[nuspec] +
                              "_" + selstr);

              c1->Clear();

              TPad *ptop = MakeRatioTopPadTopLegend();
              ptop->AppendPad();
              ptop->cd();

              NEUT->GetYaxis()->SetTitle("A.U.");
              NEUT->GetXaxis()->SetLabelSize(0);
              NEUT->GetXaxis()->SetTitleSize(0);
              NEUT->GetYaxis()->SetNdivisions(505);
              NEUT->GetXaxis()->SetNdivisions(505);
              StyleTH1Line(NEUT, kT2KRed, 2, 1);
              StyleTH1Line(NEUT_untuned, kT2KRed, 2, 2);
              StyleTH1Line(GENIE, kNOvABlue, 2, 1);
              StyleTH1Line(GENIE_untuned, kNOvABlue, 2, 2);

              DrawTH1s({NEUT, NEUT_untuned, GENIE, GENIE_untuned}, "EHIST");

              c1->cd();
              TPad *pbot = MakeRatioBottomPadTopLegend();
              pbot->AppendPad();
              pbot->cd();

              for (auto &h : {NEUT_untuned, GENIE, GENIE_untuned}) {
                h->Divide(NEUT);
              }
              NEUT_untuned->GetYaxis()->SetTitle("Ratio to NEUT+BANFF");
              NEUT_untuned->GetYaxis()->SetLabelSize(0.06);
              NEUT_untuned->GetYaxis()->SetTitleSize(0.06);
              NEUT_untuned->GetXaxis()->SetLabelSize(0.06);
              NEUT_untuned->GetXaxis()->SetTitleSize(0.06);
              NEUT_untuned->GetYaxis()->SetNdivisions(505);
              NEUT_untuned->GetXaxis()->SetNdivisions(505);

              StyleTH1Line(Inputs1D, kMSUGreen, 3, 5);

              DrawTH1s({Inputs1D, NEUT_untuned, GENIE, GENIE_untuned}, "HIST");

              c1->cd();
              auto *leg = MakeTopLegend();
              leg->SetTextSize(0.04);

              leg->AddEntry(NEUT, "NEUT+BANFF", "l");
              leg->AddEntry(NEUT_untuned, "NEUT", "l");
              leg->AddEntry(GENIE, "GENIE+NOvA2020", "l");
              leg->AddEntry(GENIE_untuned, "GENIE", "l");
              leg->AddEntry(Inputs1D, "1D Input Ratio", "l");

              leg->Draw();

              TLatex ltx;
              ltx.SetTextAlign(31);
              ltx.DrawLatexNDC(0.9, 0.75,
                               (std::string(detstr) + " " + selstr + " " +
                                projstr + " " + tgtstr + " " +
                                t2knova::all_nuspecies[nuspec])
                                   .c_str());

              c1->Print("T2KNOvATunePreds.pdf");
            }
          }
        }
        if (selstr == "totxsecs") { // don't duplicate for multiple projections
                                    // for the totxsec plot
          break;
        }
      }
    }
  }
  c1->Print("T2KNOvATunePreds.pdf]");
}