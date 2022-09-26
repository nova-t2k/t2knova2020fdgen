#include "../T2KNOvA/ROOTHelper.hxx"
#include "../include/colordef.h"
#include "../include/plotutils.h"

#include "../T2KNOvA/TrueSelectionHelper.hxx"

#include "TLatex.h"

void plot() {

  DeclareColors();

  std::vector<t2knova::selection> selections = {
      t2knova::kCCInc,  t2knova::kCCInc_RW,  t2knova::kCC0pi,  t2knova::kCC1cpi,
      t2knova::kCC1pi0, t2knova::kCCmultipi, t2knova::kCCOther};

  std::string det = "NOvAND";
  // det = "ND280";
  std::string gen = "GENIE";
  // gen = "NEUT";

  std::string tune = "2020";
  // tune = "BANFF_PRE";
  tune = "Generated";

  std::vector<std::unique_ptr<TH1>> Valids;
  std::vector<std::unique_ptr<TH1>> Hists;

  auto canv = MakeCanvasTopLegend();
  canv->Print("CCIncSums.pdf[");

  auto leg = MakeTopLegend();
  leg->SetTextSize(0.04);

  std::string valid_file = "FakeDataValid.root";
  // valid_file = "FakeDataHists_1.root";
  std::string valid_hstub = "Enu_";
  std::string hists_file = "FakeDataHists.root";
  // hists_file = "FakeDataHists_1.root";
  std::string hists_hstub = "Enu_";

  int cwheel = 0;
  for (auto s : selections) {
    std::string vname = det + "/" + gen + "/" + tune + "/CH/numu/" +
                        valid_hstub + t2knova::SelectionList[s];
    Valids.push_back(GetTH1(valid_file, vname));

    std::string hname = det + "/" + gen + "/" + tune + "/CH/numu/" +
                        hists_hstub + t2knova::SelectionList[s];
    Hists.push_back(GetTH1(hists_file, hname));

    if (!Valids.back()) {
      std::cout << "[ERROR]: Failed to read "
                << det + "/" + gen + "/" + tune + "/CH/numu/Enu_" +
                       t2knova::SelectionList[s]
                << " from " << valid_file << std::endl;
      abort();
    }

    if (!Hists.back()) {
      std::cout << "[ERROR]: Failed to read "
                << det + "/" + gen + "/" + tune + "/CH/numu/Enu_" +
                       t2knova::SelectionList[s]
                << " from " << hists_file << std::endl;
      abort();
    }

    Hists.back()->GetXaxis()->SetRangeUser(0, 5);
    Hists.back()->SetTitle("");

    Valids.back()->GetXaxis()->SetRangeUser(0, 5);
    Valids.back()->SetTitle("");

    Valids.back()->SetLineWidth(2);
    Valids.back()->SetLineColor(SORNMutedWheel[cwheel]);
    Hists.back()->SetLineWidth(2);
    Hists.back()->SetLineColor(SORNMutedWheel[cwheel]);

    leg->AddEntry(Hists.back().get(), t2knova::SelectionList[s].c_str(), "l");

    cwheel++;
  }

  auto pt = MakeRatioTopPadTopLegend();
  pt->AppendPad();
  auto pb = MakeRatioBottomPadTopLegend();
  pb->AppendPad();

  std::vector<std::unique_ptr<TH1>> ValidClones;
  {

    pt->cd();
    DrawTH1s(Hists, "HIST");

    for (int i = 0; i < Hists.size(); ++i) {
      ValidClones.push_back(Clone(Valids[i]));
      ValidClones[i]->SetLineStyle(2);
    }
    DrawTH1s(ValidClones, "HIST", false);

    for (int i = 0; i < Hists.size(); ++i) {
      ValidClones[i]->Divide(Hists[i].get());
      ValidClones[i]->SetLineStyle(1);
    }

    pb->cd();
    ValidClones.front()->GetYaxis()->SetRangeUser(0.95, 1.05);
    ValidClones.front()->GetYaxis()->SetTitle("Valids/Hists");
    DrawTH1s(ValidClones, "HIST", true, false);

    canv->cd();

    leg->Draw();

    TLatex ltx;
    ltx.SetTextSize(0.05);
    ltx.DrawLatexNDC(0.35, 0.825, (gen + " " + det + " " + tune).c_str());

    canv->Print("CCIncSums.pdf");
  }

  {

    pt->Clear();
    pt->cd();
    for (int i = (Hists.size() - 2); i >= 2; --i) {
      Hists[i]->Add(Hists[i + 1].get());
      Valids[i]->Add(Valids[i + 1].get());
    }

    Hists[2]->SetLineStyle(2);
    Valids[2]->SetLineStyle(2);

    DrawTH1s(Hists, "HIST");

    pb->Clear();
    pb->cd();
    DrawTH1s(Valids, "HIST");

    canv->cd();

    TLatex ltx;
    ltx.SetTextSize(0.05);
    ltx.DrawLatexNDC(0.35, 0.825, (gen + " " + det + " " + tune).c_str());
    ltx.DrawLatexNDC(0.7, 0.725, "Hists");
    ltx.DrawLatexNDC(0.7, 0.25, "Valids");

    canv->Print("CCIncSums.pdf");
  }

  canv->Print("CCIncSums.pdf]");
}