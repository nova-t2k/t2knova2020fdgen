#include "T2KNOvAFakeDataHelper.hxx"
#include "T2KNOvATrueSelectionHelper.hxx"
#include "TLatex.h"
#include "colordef.h"
#include "plotutils.h"

int main(int argc, char const *argv[]) {
  DeclareColors();

  std::unique_ptr<TFile> F_FDH(new TFile(argv[1], "READ"));

  TCanvas *c1 = MakeCanvasTopLegend();
  c1->SetName("T2KNOvATunePreds");
  c1->Print("T2KNOvATunePreds.pdf[");

  std::vector<std::string> selstrs = t2knova::SelectionList;
  selstrs.push_back("totxsecs");

  bool domodes = true;
  bool do0pi = true;

  TCanvas *c2 = nullptr;
  if (domodes) {
    c2 = MakeCanvasTopLegend();
    c2->SetName("T2KNOvATunePreds_modes");

    c2->Print("T2KNOvATunePreds_modes.pdf[");
  }

  for (auto detstr : {"ND280"}) {
    for (auto selstr : selstrs) {
      for (auto tgtstr : {"C", "O", "H"}) {
        for (auto nuspec :
             {t2knova::kNuMu, t2knova::kNuMub, t2knova::kNuE, t2knova::kNuEb}) {
          for (auto _projstr : {
                   "Enu",
                   "Q2",
                   "EGamma",
                   "Enu_Frac",
               }) {
            std::cout << detstr << ":" << selstr << ":" << tgtstr << ":"
                      << t2knova::all_nuspecies[nuspec] << ":" << _projstr
                      << std::endl;

            for (int mode = 0; mode < (domodes ? 120 : 1); ++mode) {
              int lmode = mode - (domodes ? 60 : 0);
              if (selstr == "totxsecs") {
                std::string mode_str = "";
                std::string projstr = _projstr;

                if (lmode != 0) {
                  mode_str = "_Mode_";
                  mode_str = mode_str + (lmode < 0 ? "m" : "") +
                             std::to_string(std::abs(lmode));
                }

                std::unique_ptr<TH1> NEUT =
                    GetTH1(F_FDH, std::string("NEUT/") + detstr + "/" + tgtstr +
                                      "/" + t2knova::all_nuspecies[nuspec] +
                                      "/" + selstr + mode_str);

                if (NEUT) {
                  std::cout << "Read: "
                            << std::string("NEUT/") + detstr + "/" + tgtstr +
                                   "/" + t2knova::all_nuspecies[nuspec] + "/" +
                                   selstr + mode_str
                            << std::endl;
                }

                std::unique_ptr<TH1> NEUT_untuned =
                    GetTH1(F_FDH, std::string("NEUT/") + detstr + "/" + tgtstr +
                                      "/" + t2knova::all_nuspecies[nuspec] +
                                      "/" + selstr + "_untuned" + mode_str);

                std::unique_ptr<TH1> GENIE = GetTH1(
                    F_FDH, std::string("GENIE/") + detstr + "/" + tgtstr + "/" +
                               t2knova::all_nuspecies[nuspec] + "/" + selstr +
                               mode_str);

                std::unique_ptr<TH1> GENIE_untuned = GetTH1(
                    F_FDH, std::string("GENIE/") + detstr + "/" + tgtstr + "/" +
                               t2knova::all_nuspecies[nuspec] + "/" + selstr +
                               "_untuned" + mode_str);

                if (!NEUT || !NEUT_untuned || !GENIE || !GENIE_untuned) {
                  continue;
                }

                if ((IntegralTH(NEUT) + IntegralTH(NEUT_untuned) +
                     IntegralTH(GENIE) + IntegralTH(GENIE_untuned)) == 0) {
                  continue;
                }

                std::unique_ptr<TH1> &first = NEUT ? NEUT : GENIE;

                c1->cd();

                first->GetYaxis()->SetTitle("#sigma_{tot.}");
                StyleTH1Line(NEUT, kT2KRed, 2, 1);
                StyleTH1Line(NEUT_untuned, kT2KRed, 2, 2);
                StyleTH1Line(GENIE, kNOvABlue, 2, 1);
                StyleTH1Line(GENIE_untuned, kNOvABlue, 2, 2);

                auto *leg = MakeTopLegend();
                leg->SetTextSize(0.04);

                if (NEUT) {
                  leg->AddEntry(NEUT.get(), "NEUT+BANFF", "l");
                  leg->AddEntry(NEUT_untuned.get(), "NEUT", "l");
                }
                if (GENIE) {
                  leg->AddEntry(GENIE.get(), "GENIE+NOvA2020", "l");
                  leg->AddEntry(GENIE_untuned.get(), "GENIE", "l");
                }

                DrawTH1s(
                    std::vector<std::reference_wrapper<std::unique_ptr<TH1>>>{
                        NEUT, NEUT_untuned, GENIE, GENIE_untuned},
                    "EHIST");
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
                std::string projstr = _projstr;

                if ((_projstr == "Enu_Frac")) {
                  if (((selstr == "CCInc") || (selstr == "NCInc"))) {
                    continue;
                  }
                  projstr = "Enu";
                }

                std::string incsel = (selstr[0] == 'C') ? "CCInc" : "NCInc";

                if (lmode != 0) {
                  mode_str = "_Mode_";
                  mode_str = mode_str + (lmode < 0 ? "m" : "") +
                             std::to_string(std::abs(lmode));
                }

                std::unique_ptr<TH1> NEUT =
                    GetTH1(F_FDH, std::string("NEUT/") + detstr + "/" + tgtstr +
                                      "/" + t2knova::all_nuspecies[nuspec] +
                                      "/" + projstr + "_" + selstr + mode_str);

                if (NEUT) {
                  std::cout << "Read: "
                            << std::string("NEUT/") + detstr + "/" + tgtstr +
                                   "/" + t2knova::all_nuspecies[nuspec] + "/" +
                                   projstr + "_" + selstr + mode_str
                            << std::endl;
                } else {
                  continue;
                }

                std::unique_ptr<TH1> NEUT_untuned = GetTH1(
                    F_FDH, std::string("NEUT/") + detstr + "/" + tgtstr + "/" +
                               t2knova::all_nuspecies[nuspec] + "/" + projstr +
                               "_untuned_" + selstr + mode_str);

                std::unique_ptr<TH1> GENIE = GetTH1(
                    F_FDH, std::string("GENIE/") + detstr + "/" + tgtstr + "/" +
                               t2knova::all_nuspecies[nuspec] + "/" + projstr +
                               "_" + selstr + mode_str);

                std::unique_ptr<TH1> GENIE_untuned = GetTH1(
                    F_FDH, std::string("GENIE/") + detstr + "/" + tgtstr + "/" +
                               t2knova::all_nuspecies[nuspec] + "/" + projstr +
                               "_untuned_" + selstr + mode_str);

                if (!NEUT || !NEUT_untuned || !GENIE || !GENIE_untuned) {
                  continue;
                }

                if ((_projstr == "Enu_Frac")) {
                  std::unique_ptr<TH1> NEUT_incsel = GetTH1(
                      F_FDH, std::string("NEUT/") + detstr + "/" + tgtstr +
                                 "/" + t2knova::all_nuspecies[nuspec] + "/" +
                                 projstr + "_" + incsel + mode_str);
                  std::unique_ptr<TH1> NEUT_untuned_incsel = GetTH1(
                      F_FDH, std::string("NEUT/") + detstr + "/" + tgtstr +
                                 "/" + t2knova::all_nuspecies[nuspec] + "/" +
                                 projstr + "_untuned_" + incsel + mode_str);

                  std::unique_ptr<TH1> GENIE_incsel = GetTH1(
                      F_FDH, std::string("GENIE/") + detstr + "/" + tgtstr +
                                 "/" + t2knova::all_nuspecies[nuspec] + "/" +
                                 projstr + "_" + incsel + mode_str);

                  std::unique_ptr<TH1> GENIE_untuned_incsel = GetTH1(
                      F_FDH, std::string("GENIE/") + detstr + "/" + tgtstr +
                                 "/" + t2knova::all_nuspecies[nuspec] + "/" +
                                 projstr + "_untuned_" + incsel + mode_str);

                  if (!NEUT_incsel || !NEUT_untuned_incsel || !GENIE_incsel ||
                      !GENIE_untuned_incsel) {
                    continue;
                  }

                  NEUT->Divide(NEUT_incsel.get());
                  NEUT_untuned->Divide(NEUT_untuned_incsel.get());
                  GENIE->Divide(GENIE_incsel.get());
                  GENIE_untuned->Divide(GENIE_untuned_incsel.get());

                  NEUT->Scale(100);
                  NEUT_untuned->Scale(100);
                  GENIE->Scale(100);
                  GENIE_untuned->Scale(100);
                  for (int x = 0; x < NEUT->GetXaxis()->GetNbins(); ++x) {
                    for (int y = 0; y < NEUT->GetYaxis()->GetNbins(); ++y) {
                      NEUT->SetBinError(x + 1, y + 1, 0);
                      NEUT_untuned->SetBinError(x + 1, y + 1, 0);
                      GENIE->SetBinError(x + 1, y + 1, 0);
                      GENIE_untuned->SetBinError(x + 1, y + 1, 0);
                    }
                  }
                }

                if ((IntegralTH(NEUT) + IntegralTH(NEUT_untuned) +
                     IntegralTH(GENIE) + IntegralTH(GENIE_untuned)) == 0) {
                  continue;
                }

                if (domodes && (lmode != 0)) {
                  c2->cd();
                  c2->Clear();
                } else {
                  c1->cd();
                  c1->Clear();
                }

                TPad *ptop = MakeRatioTopPadTopLegend();
                ptop->AppendPad();
                ptop->cd();

                std::unique_ptr<TH1> &first =
                    NEUT ? NEUT_untuned : GENIE_untuned;
                std::unique_ptr<TH1> &other =
                    NEUT ? GENIE_untuned : NEUT_untuned;

                if ((_projstr == "Enu_Frac")) {
                  first->GetYaxis()->SetTitle(
                      "Fractional XSec(#it{E}_{#nu}) %");
                } else {
                  first->GetYaxis()->SetTitle("A.U.");
                }
                first->SetTitle("");
                first->GetXaxis()->SetLabelSize(0);
                first->GetXaxis()->SetTitleSize(0);
                first->GetYaxis()->SetNdivisions(505);
                first->GetXaxis()->SetNdivisions(505);
                StyleTH1Line(NEUT, kT2KRed, 2, 1);
                StyleTH1Line(NEUT_untuned, kT2KRed, 2, 2);
                StyleTH1Line(GENIE, kNOvABlue, 2, 1);
                StyleTH1Line(GENIE_untuned, kNOvABlue, 2, 2);

                DrawTH1s(
                    std::vector<std::reference_wrapper<std::unique_ptr<TH1>>>{
                        first, other, NEUT, GENIE},
                    "EHIST");

                if (domodes && (lmode != 0)) {
                  c2->cd();
                } else {
                  c1->cd();
                }
                TPad *pbot = MakeRatioBottomPadTopLegend();
                pbot->AppendPad();
                pbot->cd();

                for (auto &h :
                     std::vector<std::reference_wrapper<std::unique_ptr<TH1>>>{
                         other, NEUT, GENIE}) {
                  if (!h.get()) {
                    continue;
                  }
                  h.get()->Divide(first.get());
                }
                other->GetYaxis()->SetTitle("Ratio to NEUT");
                if (!NEUT) {
                  other->GetYaxis()->SetTitle("Ratio to GENIE");
                }
                other->SetTitle("");
                other->GetYaxis()->SetLabelSize(0.06);
                other->GetYaxis()->SetTitleSize(0.06);
                other->GetXaxis()->SetLabelSize(0.06);
                other->GetXaxis()->SetTitleSize(0.06);
                other->GetYaxis()->SetNdivisions(505);
                other->GetXaxis()->SetNdivisions(505);

                DrawTH1s(
                    std::vector<std::reference_wrapper<std::unique_ptr<TH1>>>{
                        other, NEUT, GENIE},
                    "HIST");

                if (domodes && (lmode != 0)) {
                  c2->cd();
                } else {
                  c1->cd();
                }
                auto *leg = MakeTopLegend();
                leg->SetTextSize(0.04);
                if (NEUT) {
                  leg->AddEntry(NEUT.get(), "NEUT+BANFF", "l");
                  leg->AddEntry(NEUT_untuned.get(), "NEUT", "l");
                }
                if (GENIE) {
                  leg->AddEntry(GENIE.get(), "GENIE+NOvA2020", "l");
                  leg->AddEntry(GENIE_untuned.get(), "GENIE", "l");
                }

                leg->Draw();

                TLatex ltx;
                ltx.SetTextAlign(31);
                ltx.SetTextSize(0.04);
                ltx.DrawLatexNDC(
                    0.9, 0.8,
                    (std::string(detstr) + " " + selstr + " " + _projstr + " " +
                     tgtstr + " " + t2knova::all_nuspecies[nuspec])
                        .c_str());
                if (lmode != 0) {
                  std::string mode_nice_str = " Mode == ";
                  mode_nice_str = mode_nice_str + std::to_string(lmode);

                  ltx.DrawLatexNDC(0.9, 0.65, mode_nice_str.c_str());
                }

                if (domodes && (lmode != 0)) {
                  c2->Print("T2KNOvATunePreds_modes.pdf");
                } else {
                  c1->Print("T2KNOvATunePreds.pdf");
                }
              }
            }
            if (selstr == "totxsecs") {  // don't duplicate for multiple
                                         // projections for the totxsec plot
              break;
            }
          }
        }
      }
    }
  }
  c1->Print("T2KNOvATunePreds.pdf]");
  if (domodes) {
    c2->Print("T2KNOvATunePreds_modes.pdf]");
  }
}