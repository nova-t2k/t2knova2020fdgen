#include "T2KNOvAFakeDataHelper.hxx"

#include "colordef.h"
#include "plotutils.h"

#include "TLatex.h"

int main(int argc, char const *argv[]) {
  DeclareColors();

  // TFile F_FDH("FakeDataHists.root", "READ");
  TFile F_FDH(argv[1], "READ");
  // TFile F_FDI("FakeDataInputs.root", "READ");

  TCanvas *c1 = MakeCanvasTopLegend();
  c1->Print("T2KNOvATunePreds.pdf[");

  for (auto detstr : {"ND280"}) {
    for (auto selstr :
         {"CCInc", "CC0Pi", "CC1CPi", "CC1Pi0", "CCOther", "totxsecs"}) {
      for (auto projstr : {"Enu", "Q2"}) {
        for (auto tgtstr : {"C", "O", "H"}) {
          for (auto nuspec : {t2knova::kNuMu, t2knova::kNuMub, t2knova::kNuE,
                              t2knova::kNuEb}) {
            for (int mode = 0; mode < 120; ++mode) {
              int lmode = mode - 60;
              if (selstr == "totxsecs") {

                std::string mode_str = "";

                if (mode != 0) {
                  mode_str = "_Mode_";
                  mode_str = mode_str + (lmode < 0 ? "m" : "") +
                             std::to_string(std::abs(lmode));
                }

                TH1 *NEUT = GetTH1(&F_FDH, std::string("NEUT/") + detstr + "/" +
                                               tgtstr + "/" +
                                               t2knova::all_nuspecies[nuspec] +
                                               "/" + selstr + mode_str);

                TH1 *NEUT_untuned = GetTH1(
                    &F_FDH, std::string("NEUT/") + detstr + "/" + tgtstr + "/" +
                                t2knova::all_nuspecies[nuspec] + "/" + selstr +
                                "_untuned" + mode_str);

                TH1 *GENIE = GetTH1(&F_FDH, std::string("GENIE/") + detstr +
                                                "/" + tgtstr + "/" +
                                                t2knova::all_nuspecies[nuspec] +
                                                "/" + selstr + mode_str);

                TH1 *GENIE_untuned = GetTH1(
                    &F_FDH, std::string("GENIE/") + detstr + "/" + tgtstr +
                                "/" + t2knova::all_nuspecies[nuspec] + "/" +
                                selstr + "_untuned" + mode_str);

                if (!NEUT || !NEUT_untuned || !GENIE || !GENIE_untuned) {
                  continue;
                }
                
                if ((IntegralTH(NEUT) + IntegralTH(NEUT_untuned) +
                     IntegralTH(GENIE) + IntegralTH(GENIE_untuned)) == 0) {
                  continue;
                }

                TH1 *first = NEUT;
                if (!NEUT) {
                  first = GENIE;
                }

                c1->Clear();
                first->GetYaxis()->SetTitle("#sigma_{tot.}");
                StyleTH1Line(NEUT, kT2KRed, 2, 1);
                StyleTH1Line(NEUT_untuned, kT2KRed, 2, 2);
                StyleTH1Line(GENIE, kNOvABlue, 2, 1);
                StyleTH1Line(GENIE_untuned, kNOvABlue, 2, 2);

                auto *leg = MakeTopLegend();
                leg->SetTextSize(0.04);

                if (NEUT) {
                  leg->AddEntry(NEUT, "NEUT+BANFF", "l");
                  leg->AddEntry(NEUT_untuned, "NEUT", "l");
                }
                if (GENIE) {
                  leg->AddEntry(GENIE, "GENIE+NOvA2020", "l");
                  leg->AddEntry(GENIE_untuned, "GENIE", "l");
                }

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

                std::string mode_str = "";

                if (mode != 0) {
                  mode_str = "_Mode_";
                  mode_str = mode_str + (lmode < 0 ? "m" : "") +
                             std::to_string(std::abs(lmode));
                }

                TH1 *NEUT = GetTH1(
                    &F_FDH, std::string("NEUT/") + detstr + "/" + tgtstr + "/" +
                                t2knova::all_nuspecies[nuspec] + "/" + projstr +
                                "_" + selstr + mode_str);

                TH1 *NEUT_untuned = GetTH1(
                    &F_FDH, std::string("NEUT/") + detstr + "/" + tgtstr + "/" +
                                t2knova::all_nuspecies[nuspec] + "/" + projstr +
                                "_untuned_" + selstr + mode_str);

                TH1 *GENIE = GetTH1(
                    &F_FDH, std::string("GENIE/") + detstr + "/" + tgtstr +
                                "/" + t2knova::all_nuspecies[nuspec] + "/" +
                                projstr + "_" + selstr + mode_str);

                TH1 *GENIE_untuned = GetTH1(
                    &F_FDH, std::string("GENIE/") + detstr + "/" + tgtstr +
                                "/" + t2knova::all_nuspecies[nuspec] + "/" +
                                projstr + "_untuned_" + selstr + mode_str);

                if (!NEUT && !NEUT_untuned && !GENIE && !GENIE_untuned) {
                  continue;
                }

                if ((IntegralTH(NEUT) + IntegralTH(NEUT_untuned) +
                     IntegralTH(GENIE) + IntegralTH(GENIE_untuned)) == 0) {
                  continue;
                }

                // TH1 *Inputs1D = GetTH1(
                //     &F_FDI, std::string("t2knd_to_nova_") + projstr + "_" +
                //                 tgtstr + "_" + t2knova::all_nuspecies[nuspec]
                //                 +
                //                 "_" + selstr + mode_str);

                c1->Clear();

                TPad *ptop = MakeRatioTopPadTopLegend();
                ptop->AppendPad();
                ptop->cd();

                TH1 *first = NEUT;
                TH1 *first_next = NEUT_untuned;
                if (!NEUT) {
                  first = GENIE;
                  first_next = GENIE_untuned;
                }

                first->GetYaxis()->SetTitle("A.U.");
                first->GetXaxis()->SetLabelSize(0);
                first->GetXaxis()->SetTitleSize(0);
                first->GetYaxis()->SetNdivisions(505);
                first->GetXaxis()->SetNdivisions(505);
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
                  if (!h) {
                    continue;
                  }
                  h->Divide(first);
                }
                first_next->GetYaxis()->SetTitle("Ratio to NEUT+BANFF");
                if (!NEUT) {
                  first_next->GetYaxis()->SetTitle("Ratio to GENIE+NOvA2020");
                }
                first_next->GetYaxis()->SetLabelSize(0.06);
                first_next->GetYaxis()->SetTitleSize(0.06);
                first_next->GetXaxis()->SetLabelSize(0.06);
                first_next->GetXaxis()->SetTitleSize(0.06);
                first_next->GetYaxis()->SetNdivisions(505);
                first_next->GetXaxis()->SetNdivisions(505);

                // StyleTH1Line(Inputs1D, kMSUGreen, 3, 5);

                // DrawTH1s({Inputs1D, NEUT_untuned, GENIE, GENIE_untuned},
                //          "HIST");
                DrawTH1s({NEUT_untuned, GENIE, GENIE_untuned}, "HIST");

                c1->cd();
                auto *leg = MakeTopLegend();
                leg->SetTextSize(0.04);
                if (NEUT) {
                  leg->AddEntry(NEUT, "NEUT+BANFF", "l");
                  leg->AddEntry(NEUT_untuned, "NEUT", "l");
                }
                if (GENIE) {
                  leg->AddEntry(GENIE, "GENIE+NOvA2020", "l");
                  leg->AddEntry(GENIE_untuned, "GENIE", "l");
                }
                // leg->AddEntry(Inputs1D, "1D Input Ratio", "l");

                leg->Draw();

                TLatex ltx;
                ltx.SetTextAlign(31);
                ltx.DrawLatexNDC(0.9, 0.75,
                                 (std::string(detstr) + " " + selstr + " " +
                                  projstr + " " + tgtstr + " " +
                                  t2knova::all_nuspecies[nuspec])
                                     .c_str());
                if (mode != 0) {
                  std::string mode_nice_str = " Mode == ";
                  mode_nice_str = mode_nice_str + std::to_string(lmode);

                  ltx.DrawLatexNDC(0.9, 0.65, mode_nice_str.c_str());
                }

                c1->Print("T2KNOvATunePreds.pdf");
              }
            }
          }
        }
        if (selstr == "totxsecs") { // don't duplicate for multiple
                                    // projections for the totxsec plot
          break;
        }
      }
    }
  }
  c1->Print("T2KNOvATunePreds.pdf]");
}