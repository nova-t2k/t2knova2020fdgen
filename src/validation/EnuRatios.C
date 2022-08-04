#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"

#include <iostream>

#define CHECKH(a)                                                              \
  if (!a) {                                                                    \
    std::cout << #a << " was bad." << std::endl;                               \
    return 1;                                                                  \
  }                                                                            \
  a->SetName(#a);

#define DRAW(a, b, c)                                                          \
  a->SetLineColorAlpha(b, 0.5);                                                \
  a->SetLineWidth(2);                                                          \
  a->SetMarkerColorAlpha(b, 0);                                                \
  a->SetTitle(#a);                                                             \
  a->Draw(c);

#define DRAWS(a, b, d, c)                                                      \
  a->SetLineColor(b);                                                          \
  a->SetLineWidth(2);                                                          \
  a->SetLineStyle(d);                                                          \
  a->SetMarkerColorAlpha(b, 0);                                                \
  a->SetTitle(#a);                                                             \
  a->Draw(c);

int main() {

  TFile FDHists("FDSInputs/FakeDataHists.root");
  TFile FDInputs_Generated("FDSInputs/FakeDataInputs_FromGenerated.root");
  TFile FDInputs_Tuned("FDSInputs/FakeDataInputs_FromTuned.root");

  TFile FDValidationHists_Generated(
      "FDSValid_FromGenerated/FakeDataValid.root");
  TFile FDValidationHists_Tuned("FDSValid_FromTuned/FakeDataValid.root");

  TH1 *FDHist_Enu_Generated =
      FDHists.Get<TH1>("NEUT/ND280/C/numu/Generated/Enu_CC0pi");
  CHECKH(FDHist_Enu_Generated);
  TH1 *FDHist_Enu_BANFF_POST =
      FDHists.Get<TH1>("NEUT/ND280/C/numu/BANFF_POST/Enu_CC0pi");
  CHECKH(FDHist_Enu_BANFF_POST);
  TH1 *FDHist_Enu_2020 = FDHists.Get<TH1>("GENIE/ND280/C/numu/2020/Enu_CC0pi");
  CHECKH(FDHist_Enu_2020);

  TH1 *FDRatio_Enu_FromGenerated =
      FDInputs_Generated.Get<TH1>("Generated_to_2020/Enu/C/numu/CC0pi");
  CHECKH(FDRatio_Enu_FromGenerated);
  TH1 *FDRatio_Enu_FromTuned =
      FDInputs_Tuned.Get<TH1>("BANFF_POST_to_2020/Enu/C/numu/CC0pi");
  CHECKH(FDRatio_Enu_FromTuned);

  TH1 *FDValidHist_Enu_2020_Generated =
      FDValidationHists_Generated.Get<TH1>("ND280/GENIE/2020/C/numu/Enu_CC0pi");
  CHECKH(FDValidHist_Enu_2020_Generated);
  TH1 *FDValidHist_Enu_2020_Tuned =
      FDValidationHists_Tuned.Get<TH1>("ND280/GENIE/2020/C/numu/Enu_CC0pi");
  CHECKH(FDValidHist_Enu_2020_Tuned);

  TH1 *FDValidHist_Enu_Generated_Generated =
      FDValidationHists_Generated.Get<TH1>(
          "ND280/NEUT/Generated/C/numu/Enu_CC0pi");
  CHECKH(FDValidHist_Enu_Generated_Generated);
  TH1 *FDValidHist_Enu_BANFF_POST_Tuned = FDValidationHists_Tuned.Get<TH1>(
      "ND280/NEUT/BANFF_POST/C/numu/Enu_CC0pi");
  CHECKH(FDValidHist_Enu_BANFF_POST_Tuned);

  TH1 *FDValidHist_Enu_RW2020_Generated = FDValidationHists_Generated.Get<TH1>(
      "ND280/NEUT/ReWeighted_to_2020/C/numu/Enu_CC0pi");
  CHECKH(FDValidHist_Enu_RW2020_Generated);
  TH1 *FDValidHist_Enu_RW2020_Tuned = FDValidationHists_Tuned.Get<TH1>(
      "ND280/NEUT/ReWeighted_to_2020/C/numu/Enu_CC0pi");
  CHECKH(FDValidHist_Enu_RW2020_Tuned);

  TCanvas c1("c1", "");
  c1.Print("EnuCheck.pdf[");

  FDHist_Enu_2020->GetXaxis()->SetRangeUser(0, 3.5);
  DRAW(FDHist_Enu_2020, kBlue, "HIST");
  DRAW(FDHist_Enu_Generated, kBlack, "SAMEHIST");
  DRAW(FDHist_Enu_BANFF_POST, kRed, "SAMEHIST");

  c1.BuildLegend(0.3, 0.6, 0.9, 0.9);
  c1.Print("EnuCheck.pdf");
  FDHist_Enu_Generated->Divide(FDHist_Enu_2020, FDHist_Enu_Generated);
  FDHist_Enu_BANFF_POST->Divide(FDHist_Enu_2020, FDHist_Enu_BANFF_POST);

  FDHist_Enu_Generated->GetYaxis()->SetRangeUser(0, 3);
  FDHist_Enu_Generated->GetXaxis()->SetRangeUser(0, 3.5);
  DRAW(FDHist_Enu_Generated, kBlack, "HIST");
  DRAW(FDHist_Enu_BANFF_POST, kRed, "SAMEHIST");
  c1.BuildLegend(0.3, 0.6, 0.9, 0.9);
  c1.Print("EnuCheck.pdf");

  DRAW(FDHist_Enu_Generated, kBlack, "HIST");
  DRAW(FDHist_Enu_BANFF_POST, kRed, "SAMEHIST");
  DRAWS(FDRatio_Enu_FromGenerated, kBlack, 2, "SAMEHIST");
  DRAWS(FDRatio_Enu_FromTuned, kRed, 2, "SAMEHIST");
  c1.BuildLegend(0.3, 0.6, 0.9, 0.9);
  c1.Print("EnuCheck.pdf");

  DRAW(FDHist_Enu_2020, kBlue, "HIST");
  DRAW(FDValidHist_Enu_2020_Generated, kGreen, "SAMEHIST");
  DRAWS(FDValidHist_Enu_2020_Tuned, kGreen, 2, "SAMEHIST");
  c1.BuildLegend(0.3, 0.6, 0.9, 0.9);
  c1.Print("EnuCheck.pdf");

  FDHist_Enu_Generated =
      FDHists.Get<TH1>("NEUT/ND280/C/numu/Generated/Enu_CC0pi");
  CHECKH(FDHist_Enu_Generated);
  FDHist_Enu_BANFF_POST =
      FDHists.Get<TH1>("NEUT/ND280/C/numu/BANFF_POST/Enu_CC0pi");
  CHECKH(FDHist_Enu_BANFF_POST);

  FDHist_Enu_BANFF_POST->GetXaxis()->SetRangeUser(0, 3.5);
  DRAW(FDHist_Enu_BANFF_POST, kGreen, "HIST");
  DRAW(FDHist_Enu_Generated, kBlue, "SAMEHIST");
  DRAWS(FDValidHist_Enu_BANFF_POST_Tuned, kGreen, 2, "SAMEHIST");
  DRAWS(FDValidHist_Enu_Generated_Generated, kBlue, 2, "SAMEHIST");
  c1.BuildLegend(0.3, 0.6, 0.9, 0.9);
  c1.Print("EnuCheck.pdf");

  DRAW(FDHist_Enu_2020, kBlue, "HIST");
  DRAW(FDValidHist_Enu_RW2020_Generated, kMagenta, "SAMEHIST");
  DRAWS(FDValidHist_Enu_RW2020_Tuned, kMagenta, 2, "SAMEHIST");
  c1.BuildLegend(0.3, 0.6, 0.9, 0.9);
  c1.Print("EnuCheck.pdf");

  FDHist_Enu_Generated->GetXaxis()->SetRangeUser(0, 3.5);
  FDValidHist_Enu_RW2020_Generated->Divide(FDValidHist_Enu_RW2020_Generated,
                                           FDValidHist_Enu_Generated_Generated);
  FDValidHist_Enu_RW2020_Tuned->Divide(FDValidHist_Enu_RW2020_Tuned,
                                       FDValidHist_Enu_BANFF_POST_Tuned);
  FDHist_Enu_Generated->Divide(FDHist_Enu_2020, FDHist_Enu_Generated);
  FDHist_Enu_BANFF_POST->Divide(FDHist_Enu_2020, FDHist_Enu_BANFF_POST);
  DRAW(FDHist_Enu_Generated, kBlack, "HIST");
  DRAW(FDHist_Enu_BANFF_POST, kRed, "SAMEHIST");
  DRAW(FDValidHist_Enu_RW2020_Generated, kMagenta, "SAMEHIST");
  DRAWS(FDValidHist_Enu_RW2020_Tuned, kMagenta, 2, "SAMEHIST");
  c1.BuildLegend(0.3, 0.6, 0.9, 0.9);
  c1.Print("EnuCheck.pdf");

  c1.Print("EnuCheck.pdf]");
}