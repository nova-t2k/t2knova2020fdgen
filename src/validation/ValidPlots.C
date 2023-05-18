#include "T2KNOvA/FakeDataHelper.hxx"
#include "TLatex.h"
#include "colordef.h"
#include "plotutils.h"

#include "TExec.h"

using namespace t2knova;

enum Detector { kT2K, kNOvAND };

enum class ReWeightDenominator {
  kGenerated,
  kTuned,
  kTuned_BANFFPre,
  kTuned_BANFFPost
};

ReWeightDenominator denom = ReWeightDenominator::kGenerated;

enum class TargetTuneType {
  kGenerated,
  kTuned,
  kTuned_BANFFPre,
  kTuned_BANFFPost,
  kMnv1Pi,
  kNonQE
};

TargetTuneType tgttune = TargetTuneType::kTuned;

bool bymode = false;

std::vector<std::string> output_files_opened;

void OpenPDF(std::string const &fname) {
  if (std::find(output_files_opened.begin(), output_files_opened.end(),
                fname) == output_files_opened.end()) {
    TCanvas oc(fname.c_str(), "", 800, 800);
    oc.Print((fname + ".pdf[").c_str());
    output_files_opened.push_back(fname);
  }
}

void ClosePDFs() {
  for (auto &fname : output_files_opened) {
    TCanvas fc("fc", "", 800, 800);
    fc.Print((fname + ".pdf]").c_str());
  }
}

struct hblob {
  std::unique_ptr<TH1> From;
  std::unique_ptr<TH1> Target;
  std::unique_ptr<TH1> Outlier_high;
  std::unique_ptr<TH1> Outlier_low;
  std::vector<std::unique_ptr<TH1>> ReWeights;

  bool is2D;
  Detector _det;

  std::string varname;
  selection sel;

  void Load(std::unique_ptr<TFile> &fin, Detector det, std::string varname,
            nuspecies nuspec, std::string tgtstr, selection sel, int mode = 0) {
    _det = det;
    this->varname = varname;
    this->sel = sel;

    ReWeights.clear();

    std::string mode_str = "";
    if (mode != 0) {
      mode_str =
          (mode < 0 ? "_Mode_m" : "_Mode_") + std::to_string(std::abs(mode));
    }

    std::string sel_str = "";
    if (sel != t2knova::kNoPrimarySel) {
      sel_str = std::string("_") + SelectionList[sel];
    }

    std::string fromhist_name = "";
    std::string targethist_name = "";
    std::string RWhist_name = "";

    std::string tgtspecvar_str =
        "/" + tgtstr + "/" + all_nuspecies[nuspec] + "/" + varname;
    std::string selmode_str = sel_str + mode_str;

    if (det == kT2K) {

      switch (denom) {
      case ReWeightDenominator::kGenerated: {
        fromhist_name = "ND280/NEUT/Generated";
        break;
      }
      case ReWeightDenominator::kTuned_BANFFPre: {
        fromhist_name = "ND280/NEUT/BANFF_PRE";
        break;
      }
      case ReWeightDenominator::kTuned:
      case ReWeightDenominator::kTuned_BANFFPost: {
        fromhist_name = "ND280/NEUT/BANFF_POST";
        break;
      }
      }

      switch (tgttune) {
      case TargetTuneType::kGenerated: {
        targethist_name = "ND280/GENIE/Generated";
        RWhist_name = "ND280/NEUT/ReWeighted_to_Generated";
        break;
      }
      case TargetTuneType::kTuned: {
        targethist_name = "ND280/GENIE/2020";
        RWhist_name = "ND280/NEUT/ReWeighted_to_2020";
        break;
      }
      case TargetTuneType::kTuned_BANFFPre: {
        targethist_name = "ND280/NEUT/BANFF_PRE";
        RWhist_name = "ND280/NEUT/ReWeighted_to_BANFF_PRE";
        break;
      }
      case TargetTuneType::kTuned_BANFFPost: {
        targethist_name = "ND280/NEUT/BANFF_POST";
        RWhist_name = "ND280/NEUT/ReWeighted_to_BANFF_POST";
        break;
      }
      case TargetTuneType::kMnv1Pi: {
        targethist_name = "ND280/NEUT/Mnv1Pi";
        RWhist_name = "ND280/NEUT/ReWeighted_to_Mnv1Pi";
        break;
      }
      case TargetTuneType::kNonQE: {
        targethist_name = "ND280/NEUT/NonQE";
        RWhist_name = "ND280/NEUT/ReWeighted_to_NonQE";
        break;
      }
      }

    } else {

      switch (denom) {
      case ReWeightDenominator::kGenerated: {
        fromhist_name = "NOvAND/NEUT/Generated";
        break;
      }
      case ReWeightDenominator::kTuned: {
        fromhist_name = "NOvAND/GENIE/2020";
        break;
      }
      case ReWeightDenominator::kTuned_BANFFPre:
      case ReWeightDenominator::kTuned_BANFFPost: {
        std::cout << "[ERROR]: Invalid reweight denominator for NOvAND "
                     "validation plots."
                  << std::endl;
        throw;
      }
      }

      switch (tgttune) {
      case TargetTuneType::kGenerated: {
        targethist_name = "NOvAND/NEUT/Generated";
        RWhist_name = "NOvAND/GENIE/ReWeighted_to_Generated";
        break;
      }
      case TargetTuneType::kTuned_BANFFPre: {
        targethist_name = "NOvAND/NEUT/BANFF_PRE";
        RWhist_name = "NOvAND/GENIE/ReWeighted_to_BANFF_PRE";
        break;
      }
      case TargetTuneType::kTuned:
      case TargetTuneType::kTuned_BANFFPost: {
        targethist_name = "NOvAND/NEUT/BANFF_POST";
        RWhist_name = "NOvAND/GENIE/ReWeighted_to_BANFF_POST";
        break;
      }
      case TargetTuneType::kMnv1Pi: {
        targethist_name = "NOvAND/NEUT/Mnv1Pi";
        RWhist_name = "NOvAND/GENIE/ReWeighted_to_Mnv1Pi";
        break;
      }
      case TargetTuneType::kNonQE: {
        targethist_name = "NOvAND/NEUT/NonQE";
        RWhist_name = "NOvAND/GENIE/ReWeighted_to_NonQE";
        break;
      }
      }
    }

    fromhist_name += tgtspecvar_str + selmode_str;
    targethist_name += tgtspecvar_str + selmode_str;

    From = GetTH1(fin, fromhist_name, false);
    Target = GetTH1(fin, targethist_name, false);

    if (From) {
      From->SetName("FROM");
    }
    if (Target) {
      Target->SetName("TARGET");
    }

    std::vector<std::string> ReWeightHists;

    ReWeightHists.push_back(RWhist_name + tgtspecvar_str + selmode_str);

    int i = 0;
    for (auto &n : ReWeightHists) {
      ReWeights.push_back(GetTH1(fin, n, false));
      if (!ReWeights.back()) {
        std::cout << "[ERROR]: Failed to read " << n << std::endl;
        continue;
      }
      ReWeights.back()->SetName(
          (std::string("ReWeight_") + std::to_string(i++)).c_str());
      std::cout << "read " << n << " as " << ReWeights.back()->GetName()
                << std::endl;
    }

    if (ReWeights.size() && ReWeights[0]) {

      std::string Outlier_high_name =
          RWhist_name + tgtspecvar_str + "_outlier_high" + selmode_str;
      std::string Outlier_low_name =
          RWhist_name + tgtspecvar_str + "_outlier_low" + selmode_str;

      Outlier_low = GetTH1(fin, Outlier_low_name, false);
      Outlier_high = GetTH1(fin, Outlier_high_name, false);
    }
  }

  void Print1D(std::string fname, std::string title = "",
               std::string mode_line = "") {

    bool HaveBoth = bool(From) && bool(Target);
    bool HaveReWeight = ReWeights.size() && bool(ReWeights[0]);

    TH1 *First = (From ? From.get() : Target.get());

    double max_gen = GetMaximumTH1s(
        std::vector<std::reference_wrapper<std::unique_ptr<TH1>>>{From,
                                                                  Target});
    double max_rw = GetMaximumTH1s(ReWeights);
    double max = std::max(max_gen, max_rw);

    TCanvas *c1;
    TPad *p1, *p2;

    if (HaveBoth) {
      c1 = MakeCanvasTopLegend();
      p1 = MakeRatioTopPadTopLegend();
      p1->AppendPad();
      p2 = MakeRatioBottomPadTopLegend();
      p2->AppendPad();

      p1->cd();
    } else {
      c1 = MakeCanvas();
      c1->SetTopMargin(0.2);
    }

    if (From) {
      From->SetLineColor(kBlack);
      From->SetLineWidth(3);
    }

    if (Target) {
      Target->SetLineWidth(3);
      Target->SetLineColor(SORNBrightWheel[2]);
    }

    int cols[] = {SORNBrightWheel[4]};
    int lc = 0;
    for (auto &h : ReWeights) {
      if (!h) {
        lc++;
        continue;
      }
      h->SetLineWidth(4);
      h->SetLineStyle(2);
      h->SetLineColor(cols[lc++]);
    }

    First->GetYaxis()->SetRangeUser(0, max * 1.2);
    First->GetYaxis()->SetTitle("Cross Section 10^{-39}");
    First->SetTitle("");

    First->GetYaxis()->SetNdivisions(505);
    First->GetYaxis()->SetLabelSize(0.07);
    First->GetYaxis()->SetLabelOffset(0.01);
    First->GetYaxis()->SetTitleSize(0.07);
    First->GetYaxis()->SetTitleOffset(1);

    if (HaveBoth) {
      HideAxis(First->GetXaxis());
    }

    if (HaveReWeight) {
      ReWeights[0]->Draw();
    }

    if (((sel % kCCOther_QE) < kCCmultipi) &&
        ((varname == "Enu") || (varname == "ERecQE"))) {
      First->GetXaxis()->SetRangeUser(First->GetXaxis()->GetBinLowEdge(1),
                                      First->GetXaxis()->GetBinUpEdge(35));
    }
    First->DrawClone("EHIST");
    if (HaveBoth) {
      Target->DrawClone("EHISTSAME");
    }

    for (auto &h : ReWeights) {
      if (!h) {
        continue;
      }
      h->DrawClone("EHISTSAME");
    }

    if (HaveBoth) {
      p2->cd();

      From->Divide(Target.get());
      for (auto &h : ReWeights) {
        if (!h) {
          continue;
        }
        h->Divide(Target.get());
      }

      double max_rat =
          1.1 * std::max(GetMaximumTH1s(ReWeights), GetMaximumTH1s({
                                                        From,
                                                    }));
      double min_rat = std::min(GetMinimumTH1s(ReWeights), GetMinimumTH1s({
                                                               From,
                                                           }));

      Target->Divide(Target.get());

      Target->GetYaxis()->SetNdivisions(505);
      Target->GetYaxis()->SetLabelSize(0.07 * 2);
      Target->GetYaxis()->SetLabelOffset(0.01);
      Target->GetYaxis()->SetTitleSize(0.07 * 2);
      Target->GetYaxis()->SetTitleOffset(1. / 2.);

      Target->GetXaxis()->SetNdivisions(505);
      Target->GetXaxis()->SetLabelSize(0.07 * 2);
      Target->GetXaxis()->SetLabelOffset(0.01);
      Target->GetXaxis()->SetTitleSize(0.07 * 2);
      Target->GetXaxis()->SetTitleOffset(1.1);

      Target->GetYaxis()->SetRangeUser(min_rat, max_rat);
      Target->GetYaxis()->SetTitle("RW/Gen.");
      Target->SetTitle("");
      Target->SetLineStyle(2);

      if (((sel % kCCOther_QE) < kCCmultipi) &&
          ((varname == "Enu") || (varname == "ERecQE"))) {
        Target->GetXaxis()->SetRangeUser(Target->GetXaxis()->GetBinLowEdge(1),
                                         Target->GetXaxis()->GetBinUpEdge(35));
      }

      Target->Draw("HIST");
      From->DrawClone("HISTSAME");

      for (auto &h : ReWeights) {
        if (!h) {
          continue;
        }
        h->DrawClone("HISTSAME");
      }
      Target->SetLineStyle(1);
    }

    c1->cd();

    std::string FromTitle = "";
    std::string TargetTitle = "";
    std::string ReWeightTitle = "";

    if (_det == kT2K) {

      switch (denom) {
      case ReWeightDenominator::kGenerated: {
        FromTitle = "Untuned NEUT";
        break;
      }
      case ReWeightDenominator::kTuned_BANFFPre: {
        FromTitle = "BANFF Pre";
        break;
      }
      case ReWeightDenominator::kTuned:
      case ReWeightDenominator::kTuned_BANFFPost: {
        FromTitle = "BANFF Post";
        break;
      }
      }

      switch (tgttune) {
      case TargetTuneType::kGenerated: {
        TargetTitle = "Untuned GENIE";
        ReWeightTitle = "R/W to GENIE (E_{#nu}, p_{l}, #theta_{l})";
        break;
      }
      case TargetTuneType::kTuned: {
        TargetTitle = "NOvA2020";
        ReWeightTitle = "R/W to NOvA2020 (E_{#nu}, p_{l}, #theta_{l})";
        break;
      }
      case TargetTuneType::kTuned_BANFFPre: {
        TargetTitle = "BANFF Prefit";
        ReWeightTitle = "R/W to BANFFPre (E_{#nu}, p_{l}, #theta_{l})";
        break;
      }
      case TargetTuneType::kTuned_BANFFPost: {
        TargetTitle = "BANFF Postfit";
        ReWeightTitle = "R/W to BANFFPost (E_{#nu}, p_{l}, #theta_{l})";
        break;
      }
      case TargetTuneType::kMnv1Pi: {
        TargetTitle = "T2K Mnv1Pi";
        ReWeightTitle = "R/W to Mnv1Pi (E_{#nu}, p_{l}, #theta_{l})";
        break;
      }
      case TargetTuneType::kNonQE: {
        TargetTitle = "T2K NonQE";
        ReWeightTitle = "R/W to NonQE (E_{#nu}, p_{l}, #theta_{l})";
        break;
      }
      }

    } else { // novand

      switch (denom) {
      case ReWeightDenominator::kGenerated: {
        FromTitle = "Untuned GENIE";
        break;
      }
      case ReWeightDenominator::kTuned: {
        FromTitle = "NOvA2020";
        break;
      }
      }

      switch (tgttune) {
      case TargetTuneType::kGenerated: {
        TargetTitle = "Untuned NEUT";
        ReWeightTitle = "R/W to NEUT (E_{#nu} p_{T} EVisHad)";
        break;
      }
      case TargetTuneType::kTuned_BANFFPre: {
        TargetTitle = "BANFF Prefit";
        ReWeightTitle = "R/W to BANFFPre (E_{#nu} p_{T} EVisHad)";
        break;
      }
      case TargetTuneType::kTuned:
      case TargetTuneType::kTuned_BANFFPost: {
        TargetTitle = "BANFF Postfit";
        ReWeightTitle = "R/W to BANFFPost (E_{#nu} p_{T} EVisHad)";
        break;
      }
      case TargetTuneType::kMnv1Pi: {
        TargetTitle = "T2K Mnv1Pi";
        ReWeightTitle = "R/W to Mnv1Pi (E_{#nu} p_{T} EVisHad)";
        break;
      }
      case TargetTuneType::kNonQE: {
        TargetTitle = "T2K NonQE";
        ReWeightTitle = "R/W to NonQE (E_{#nu} p_{T} EVisHad)";
        break;
      }
      }
    }

    TLegend *leg = new TLegend(0, 0.81, 0.8, 1);
    leg->SetNColumns(2);
    leg->SetTextSize(0.06);
    leg->SetBorderSize(0);
    leg->SetFillStyle(4001);

    leg->AddEntry(From.get(), FromTitle.c_str(), "l");
    leg->AddEntry(Target.get(), TargetTitle.c_str(), "l");

    if (HaveReWeight) {
      leg->AddEntry(ReWeights[0].get(), ReWeightTitle.c_str(), "l");
    }

    leg->Draw();

    TLatex ttl;
    ttl.SetTextSize(0.05);
    ttl.SetTextAlign(32);
    ttl.DrawLatexNDC(0.95, 0.76, title.c_str());
    if (mode_line.size()) {
      ttl.DrawLatexNDC(0.95, 0.71, mode_line.c_str());
    }

    OpenPDF(fname);
    c1->Print((fname + ".pdf").c_str());
  }

  void Print2D(std::string fname, std::string title = "",
               std::string mode_line = "") {

    std::string FromTitle = "";
    std::string TargetTitle = "";

    if (_det == kT2K) {

      switch (denom) {
      case ReWeightDenominator::kGenerated: {
        FromTitle = "Untuned NEUT";
        break;
      }
      case ReWeightDenominator::kTuned_BANFFPre: {
        FromTitle = "BANFF Pre";
        break;
      }
      case ReWeightDenominator::kTuned:
      case ReWeightDenominator::kTuned_BANFFPost: {
        FromTitle = "BANFF Post";
        break;
      }
      }

      switch (tgttune) {
      case TargetTuneType::kGenerated: {
        TargetTitle = "Untuned GENIE";
        break;
      }
      case TargetTuneType::kTuned: {
        TargetTitle = "NOvA2020";
        break;
      }
      case TargetTuneType::kTuned_BANFFPre: {
        TargetTitle = "BANFF Prefit";
        break;
      }
      case TargetTuneType::kTuned_BANFFPost: {
        TargetTitle = "BANFF Postfit";
        break;
      }
      case TargetTuneType::kMnv1Pi: {
        TargetTitle = "T2K Mnv1Pi";
        break;
      }
      case TargetTuneType::kNonQE: {
        TargetTitle = "T2K NonQE";
        break;
      }
      }

    } else {

      switch (denom) {
      case ReWeightDenominator::kGenerated: {
        FromTitle = "Untuned GENIE";
        break;
      }
      case ReWeightDenominator::kTuned: {
        FromTitle = "NOvA2020";
        break;
      }
      }

      switch (tgttune) {
      case TargetTuneType::kGenerated: {
        TargetTitle = "Untuned NEUT";
        break;
      }
      case TargetTuneType::kTuned_BANFFPre: {
        TargetTitle = "BANFF Prefit";
        break;
      }
      case TargetTuneType::kTuned:
      case TargetTuneType::kTuned_BANFFPost: {
        TargetTitle = "BANFF Postfit";
        break;
      }
      case TargetTuneType::kMnv1Pi: {
        TargetTitle = "T2K Mnv1Pi";
        break;
      }
      case TargetTuneType::kNonQE: {
        TargetTitle = "T2K NonQE";
        break;
      }
      }
    }

    static const long unsigned int BWRPalette =
        (long unsigned int)GetBWRPalette();

    static std::string BWR_cmd = "gStyle->SetPalette(100, (int*)" +
                                 std::to_string(BWRPalette) +
                                 ");gStyle->SetNumberContours(100);";

    static TExec *ex1 = new TExec("ex1", "gStyle->SetPalette(kBird);");
    static TExec *ex2 = new TExec("ex2", BWR_cmd.c_str());

    bool HaveBoth = bool(From) && bool(Target);
    bool HaveReWeight = ReWeights.size() && bool(ReWeights[0]);

    bool HaveOutliers = bool(Outlier_low) && bool(Outlier_high);

    TCanvas *c1 = MakeCanvas();

    TLatex ttl;
    ttl.SetTextSize(0.05);
    ttl.SetTextAlign(22);

    if (From) {
      TPad *p = new TPad("pFrom", "", 0, 0.55, 0.5, 1);
      p->SetGridx();
      p->SetGridy();
      p->AppendPad();
      p->SetBottomMargin(0.08);
      p->SetLeftMargin(0.28);
      p->SetRightMargin(0.04);
      p->SetTopMargin(0.24);
      p->cd();

      From->GetYaxis()->SetNdivisions(505);
      From->GetYaxis()->SetLabelSize(0.1);
      From->GetYaxis()->SetLabelOffset(0.01);
      From->GetYaxis()->SetTitleSize(0.1);
      From->GetYaxis()->SetTitleOffset(1.35);
      HideAxis(From->GetXaxis());
      HideAxis(From->GetZaxis());
      From->GetXaxis()->SetNdivisions(505);
      From->GetYaxis()->SetNdivisions(505);

      From->SetTitle("");
      auto CL = From->DrawClone("COL");
      ex1->Draw();
      CL->Draw("COL SAME");
      p->RedrawAxis("g");

      c1->cd();
      ttl.DrawLatexNDC(0.33, 0.92, FromTitle.c_str());
    }

    if (Target) {
      TPad *p = new TPad("pTarget", "", 0.5, 0.55, 1, 1);
      p->SetGridx();
      p->SetGridy();
      p->AppendPad();
      p->SetBottomMargin(0.08);
      p->SetLeftMargin(0.04);
      p->SetRightMargin(0.28);
      p->SetTopMargin(0.24);
      p->cd();

      HideAxis(Target->GetYaxis());
      HideAxis(Target->GetXaxis());
      Target->GetXaxis()->SetNdivisions(505);
      Target->GetYaxis()->SetNdivisions(505);

      Target->GetZaxis()->SetNdivisions(505);
      Target->GetZaxis()->SetLabelSize(0.07);
      Target->GetZaxis()->SetLabelOffset(0.01);
      Target->GetZaxis()->SetTitleSize(0.07);
      Target->GetZaxis()->SetTitleOffset(1.5);

      Target->SetTitle("");
      auto CL = Target->DrawClone("COLZ");
      ex1->Draw();
      CL->Draw("COLZ SAME");

      p->RedrawAxis("g");

      c1->cd();
      ttl.DrawLatexNDC(0.7, 0.92, TargetTitle.c_str());
    }

    if (HaveBoth) {
      auto FROMCL = (TH2 *)From->Clone("FROMCL");
      FROMCL->Divide(Target.get());

      TPad *p = new TPad("pBoth", "", 0, 0.05, 0.5, 0.55);
      p->SetGridx();
      p->SetGridy();
      p->AppendPad();
      p->SetBottomMargin(0.24);
      p->SetLeftMargin(0.28);
      p->SetRightMargin(0.04);
      p->SetTopMargin(0.08);
      p->cd();

      FROMCL->GetYaxis()->SetNdivisions(505);
      FROMCL->GetYaxis()->SetLabelSize(0.1);
      FROMCL->GetYaxis()->SetLabelOffset(0.01);
      FROMCL->GetYaxis()->SetTitleSize(0.1);
      FROMCL->GetYaxis()->SetTitleOffset(1.35);

      FROMCL->GetXaxis()->SetNdivisions(505);
      FROMCL->GetXaxis()->SetLabelSize(0.1);
      FROMCL->GetXaxis()->SetLabelOffset(0.01);
      FROMCL->GetXaxis()->SetTitleSize(0.1);
      FROMCL->GetXaxis()->SetTitleOffset(1);

      FROMCL->GetZaxis()->SetRangeUser(0, 2);

      auto CL = FROMCL->DrawClone("COL");
      ex2->Draw();
      CL->Draw("COL SAME");
      p->RedrawAxis("g");

      c1->cd();
    }

    if (HaveReWeight) {
      TPad *p = new TPad("pRW", "", 0.5, 0.05, 1, 0.55);
      p->SetGridx();
      p->SetGridy();
      p->AppendPad();
      p->SetBottomMargin(0.24);
      p->SetLeftMargin(0.04);
      p->SetRightMargin(0.28);
      p->SetTopMargin(0.08);
      p->cd();

      ReWeights[0]->Divide(Target.get());
      ReWeights[0]->SetTitle("");

      ReWeights[0]->GetZaxis()->SetNdivisions(505);
      ReWeights[0]->GetZaxis()->SetLabelSize(0.07);
      ReWeights[0]->GetZaxis()->SetLabelOffset(0.01);
      ReWeights[0]->GetZaxis()->SetTitleSize(0.07);
      ReWeights[0]->GetZaxis()->SetTitleOffset(1.5);
      ReWeights[0]->GetZaxis()->SetRangeUser(0, 2);
      ReWeights[0]->GetZaxis()->SetTitle("Evt. Rate Ratio");

      ReWeights[0]->GetXaxis()->SetNdivisions(505);
      ReWeights[0]->GetXaxis()->SetLabelSize(0.1);
      ReWeights[0]->GetXaxis()->SetLabelOffset(0.01);
      ReWeights[0]->GetXaxis()->SetTitleSize(0.1);
      ReWeights[0]->GetXaxis()->SetTitleOffset(1);

      HideAxis(ReWeights[0]->GetYaxis());
      ReWeights[0]->GetXaxis()->SetNdivisions(505);
      ReWeights[0]->GetYaxis()->SetNdivisions(505);

      auto CL = ReWeights[0]->DrawClone("COLZ");
      ex2->Draw();
      CL->Draw("COLZ SAME");
      p->RedrawAxis("g");

      c1->cd();

      ttl.DrawLatexNDC(0.33, 0.55, (TargetTitle + "/" + FromTitle).c_str());
      ttl.DrawLatexNDC(0.7, 0.55, "RW/Gen.");
    }
    OpenPDF(fname);
    c1->Print((fname + ".pdf").c_str());

    if (!HaveOutliers) {
      return;
    }

    if (From) {
      if (Outlier_low) {
        TPad *p_low = new TPad("pFrom_low", "", 0, 0.55, 0.5, 1);
        p_low->SetGridx();
        p_low->SetGridy();
        p_low->AppendPad();
        p_low->SetBottomMargin(0.08);
        p_low->SetLeftMargin(0.28);
        p_low->SetRightMargin(0.04);
        p_low->SetTopMargin(0.24);
        p_low->cd();

        From->GetYaxis()->SetNdivisions(505);
        From->GetYaxis()->SetLabelSize(0.1);
        From->GetYaxis()->SetLabelOffset(0.01);
        From->GetYaxis()->SetTitleSize(0.1);
        From->GetYaxis()->SetTitleOffset(1.35);
        HideAxis(From->GetXaxis());
        HideAxis(From->GetZaxis());
        From->GetXaxis()->SetNdivisions(505);
        From->GetYaxis()->SetNdivisions(505);

        From->SetTitle("");
        auto CL_low = (TH2 *)From->DrawClone("COLZ");
        ex1->Draw();
        CL_low->Draw("COLZ SAME");

        Outlier_low->SetLineColor(kRed);
        Outlier_low->Scale(CL_low->GetMaximum() / Outlier_low->GetMaximum());
        Outlier_low->DrawClone("BOX SAME");

        c1->cd();
        ttl.DrawLatexNDC(0.33, 0.92, (FromTitle + " Low").c_str());
      }

      if (Outlier_high) {
        TPad *p_high = new TPad("pFrom_high", "", 0, 0.05, 0.5, 0.55);
        p_high->SetGridx();
        p_high->SetGridy();
        p_high->AppendPad();
        p_high->SetBottomMargin(0.24);
        p_high->SetLeftMargin(0.28);
        p_high->SetRightMargin(0.04);
        p_high->SetTopMargin(0.08);
        p_high->cd();

        From->GetYaxis()->SetNdivisions(505);
        From->GetYaxis()->SetLabelSize(0.1);
        From->GetYaxis()->SetLabelOffset(0.01);
        From->GetYaxis()->SetTitleSize(0.1);
        From->GetYaxis()->SetTitleOffset(1.35);

        From->GetXaxis()->SetNdivisions(505);
        From->GetXaxis()->SetLabelSize(0.1);
        From->GetXaxis()->SetLabelOffset(0.01);
        From->GetXaxis()->SetTitleSize(0.1);
        From->GetXaxis()->SetTitleOffset(1);

        From->GetZaxis()->SetNdivisions(505);
        From->GetZaxis()->SetLabelSize(0.07);
        From->GetZaxis()->SetLabelOffset(0.01);
        From->GetZaxis()->SetTitleSize(0.07);
        From->GetZaxis()->SetTitleOffset(1.5);

        From->SetTitle("");
        auto CL_high = (TH2 *)From->DrawClone("COLZ");
        ex1->Draw();
        CL_high->Draw("COLZ SAME");

        Outlier_high->SetLineColor(kRed);
        Outlier_high->Scale(CL_high->GetMaximum() / Outlier_high->GetMaximum());
        Outlier_high->DrawClone("BOX SAME");

        c1->cd();
        ttl.DrawLatexNDC(0.33, 0.55, (FromTitle + " High").c_str());
      }
    }

    if (Target) {
      if (Outlier_low) {
        TPad *p_low = new TPad("pTarget_low", "", 0.5, 0.55, 1, 1);
        p_low->SetGridx();
        p_low->SetGridy();
        p_low->AppendPad();
        p_low->SetBottomMargin(0.08);
        p_low->SetLeftMargin(0.04);
        p_low->SetRightMargin(0.28);
        p_low->SetTopMargin(0.24);
        p_low->cd();

        HideAxis(Target->GetYaxis());
        HideAxis(Target->GetXaxis());
        Target->GetXaxis()->SetNdivisions(505);
        Target->GetYaxis()->SetNdivisions(505);

        Target->GetZaxis()->SetNdivisions(505);
        Target->GetZaxis()->SetLabelSize(0.07);
        Target->GetZaxis()->SetLabelOffset(0.01);
        Target->GetZaxis()->SetTitleSize(0.07);
        Target->GetZaxis()->SetTitleOffset(1.5);

        Target->SetTitle("");
        auto CL_low = (TH2 *)Target->DrawClone("COLZ");
        ex1->Draw();
        CL_low->Draw("COLZ SAME");

        Outlier_low->SetLineColor(kRed);
        Outlier_low->Scale(CL_low->GetMaximum() / Outlier_low->GetMaximum());
        Outlier_low->DrawClone("BOX SAME");

        c1->cd();
        ttl.DrawLatexNDC(0.7, 0.92, (TargetTitle + " Low").c_str());
      }

      if (Outlier_high) {
        TPad *p_high = new TPad("pTarget_high", "", 0.5, 0.05, 1, 0.55);
        p_high->SetGridx();
        p_high->SetGridy();
        p_high->AppendPad();
        p_high->SetBottomMargin(0.24);
        p_high->SetLeftMargin(0.04);
        p_high->SetRightMargin(0.28);
        p_high->SetTopMargin(0.08);
        p_high->cd();

        Target->GetZaxis()->SetNdivisions(505);
        Target->GetZaxis()->SetLabelSize(0.07);
        Target->GetZaxis()->SetLabelOffset(0.01);
        Target->GetZaxis()->SetTitleSize(0.07);
        Target->GetZaxis()->SetTitleOffset(1.5);

        Target->GetXaxis()->SetNdivisions(505);
        Target->GetXaxis()->SetLabelSize(0.1);
        Target->GetXaxis()->SetLabelOffset(0.01);
        Target->GetXaxis()->SetTitleSize(0.1);
        Target->GetXaxis()->SetTitleOffset(1);

        HideAxis(Target->GetYaxis());
        Target->GetXaxis()->SetNdivisions(505);
        Target->GetYaxis()->SetNdivisions(505);

        Target->SetTitle("");
        auto CL_high = (TH2 *)Target->DrawClone("COLZ");
        ex1->Draw();
        CL_high->Draw("COLZ SAME");

        Outlier_high->SetLineColor(kRed);
        Outlier_high->Scale(CL_high->GetMaximum() / Outlier_high->GetMaximum());
        Outlier_high->DrawClone("BOX SAME");

        c1->cd();
        ttl.DrawLatexNDC(0.7, 0.55, (TargetTitle + " High").c_str());
      }
    }

    c1->Print((fname + ".pdf").c_str());
  }

  void Print(std::string fname, std::string title = "",
             std::string mode_line = "") {

    if (!From && !Target) {
      return;
    }

    TH1 *First = (From ? From.get() : Target.get());

    if (First->GetDimension() == 1) {
      Print1D(fname, title, mode_line);
    } else if (First->GetDimension() == 2) {
      Print2D(fname, title, mode_line);
    }
  }

  static void LoadAndPrint(std::unique_ptr<TFile> &fin, Detector det,
                           std::string varname, nuspecies nuspec,
                           std::string tgtstr, selection sel, std::string fname,
                           std::string title) {
    hblob h;
    h.Load(fin, det, varname, nuspec, tgtstr, sel);
    h.Print(fname, title + " " + tgtstr + " " + all_nuspecies_latex[nuspec],
            "");

    if (bymode) {
      for (int i = -60; i < 60; ++i) {
        if (i == 0) {
          continue;
        }
        h.Load(fin, det, varname, nuspec, tgtstr, sel, i);
        h.Print(fname + "_Modes",
                title + " " + tgtstr + " " + all_nuspecies_latex[nuspec],
                "True mode:" + std::to_string(i));
      }
    }
  }
};

void ValidPlots(std::string const &finname, bool DoT2K) {
  std::unique_ptr<TFile> fin(new TFile(finname.c_str()));
  if (fin->IsZombie()) {
    std::cout << "Failed to read " << finname << std::endl;
    return;
  }

  TCanvas c1("cdummy", "", 800, 800);

  for (auto tgtstr : {std::string("CH"), std::string("H2O")}) {
    for (auto nuspec : {kNuMu, kNuMub, kNuE, kNuEb}) {

      if (DoT2K) {
        hblob::LoadAndPrint(
            fin, kT2K, "SelectionXSecs", nuspec, tgtstr, kNoPrimarySel,
            std::string("Valid_ND280_SelectionXSecs_") + all_nuspecies[nuspec],
            "SelectionXSecs");

        for (int sel :
             {kCCInc,      kCCInc_RW,    kCC0pi,           kCC0pi_QE,
              kCC0pi_2p2h, kCC0pi_Other, kCC1Gamma,        kCCDeExciteGamma,
              kCCNGamma,   kCC1cpi,      kCC1pi0,          kCCmultipi,
              kCCOther,    kCCOther_QE,  kNCInc,           kNCInc_RW,
              kNC0pi,      kNC1Gamma,    kNCDeExciteGamma, kNCNGamma,
              kNC1cpi,     kNC1pi0,      kNCmultipi,       kNCOther,
              kNCOther_QE}) {

          for (auto t2k_proj : {
                   "Enu", "ERecQE", "PLep",
                   // "ThetaLep",
                   "CosThetaLep", "PThetaLep", "Q2", "EnuQ2", "EnuERecQE",
                   "EnuERecQEBias",
                   // "q0q3_low", "q0q3_high",
                   "hmfscpip", "hmfspi0p", "ncpi", "npi0", "hmfsprotonp",
                   "hmfsneutronp", "nproton", "nneutron",
                   // "EGamma",
                   // "EGamma_DeExcite",
                   "q0", "yrec", "Enuyrec",
                   // "EvWeights"
               }) {

            hblob::LoadAndPrint(fin, kT2K, t2k_proj, nuspec, tgtstr,
                                selection(sel),
                                "Valid_ND280_" + SelectionList[sel] + "_" +
                                    all_nuspecies[nuspec],
                                SelectionList[sel]);
          }
        }
      } else {
        hblob::LoadAndPrint(
            fin, kNOvAND, "SelectionXSecs", nuspec, tgtstr, kNoPrimarySel,
            std::string("Valid_NOvAND_SelectionXSecs_") + all_nuspecies[nuspec],
            "SelectionXSecs");

        for (int sel :
             {kCCInc,      kCCInc_RW,    kCC0pi,           kCC0pi_QE,
              kCC0pi_2p2h, kCC0pi_Other, kCC1Gamma,        kCCDeExciteGamma,
              kCCNGamma,   kCC1cpi,      kCC1pi0,          kCCmultipi,
              kCCOther,    kCCOther_QE,  kNCInc,           kNCInc_RW,
              kNC0pi,      kNC1Gamma,    kNCDeExciteGamma, kNCNGamma,
              kNC1cpi,     kNC1pi0,      kNCmultipi,       kNCOther,
              kNCOther_QE}) {

          for (auto t2k_proj : {
                   "Enu", "Q2", "PLep", "PtLep", "EAvHad", "EnuEAvHad",
                   "EnuEAvHadBias",
                   // "q0q3_low",
                   //  "q0q3_high",
                   "CosThetaLep", "hmfscpip", "hmfspi0p", "ncpi", "npi0",
                   // "EGamma", "EGamma_DeExcite"
                   "yrec"
               }) {

            hblob::LoadAndPrint(fin, kNOvAND, t2k_proj, nuspec, tgtstr,
                                selection(sel),
                                "Valid_NOvAND_" + SelectionList[sel] + "_" +
                                    all_nuspecies[nuspec],
                                SelectionList[sel]);
          }
        }
      }
    }
  }

  ClosePDFs();
}

void SayUsage(char const *argv[]) {
  std::cout << "[USAGE]: " << argv[0]
            << "\n"
               "\n\t-i <inp.root>            : Input file"
               "\n\t--From <Tuned|Tuned_BANFFPre|Tuned_BANFFPost> : Weight as "
               "if Tuned/Not Tuned"
               "\n\t--To   "
               "<Tuned|Generated|Tuned|Tuned_BANFFPre|Tuned_BANFFPost|Mnv1Pi|"
               "NonQE> : Weight as if Tuned/Not Tuned"
               "\n\t--NOvAND"
            << std::endl;
}

std::string input_file = "";

bool JustDoT2K = true;

void handleOpts(int argc, char const *argv[]) {
  int opt = 1;
  while (opt < argc) {
    if ((std::string(argv[opt]) == "-?") ||
        (std::string(argv[opt]) == "--help")) {
      SayUsage(argv);
      exit(0);
    } else if (std::string(argv[opt]) == "-i") {
      input_file = argv[++opt];
    } else if (std::string(argv[opt]) == "--NOvAND") {
      JustDoT2K = false;
    } else if (std::string(argv[opt]) == "--From") {
      std::string arg = std::string(argv[++opt]);

      if (arg == "Tuned") {
        denom = ReWeightDenominator::kTuned;
      } else if (arg == "Tuned_BANFFPre") {
        denom = ReWeightDenominator::kTuned_BANFFPre;
      } else if (arg == "Tuned_BANFFPost") {
        denom = ReWeightDenominator::kTuned_BANFFPost;
      } else {
        std::cout << "[ERROR]: Invalid option passed to --From " << arg
                  << std::endl;
        abort();
      }

    } else if (std::string(argv[opt]) == "--To") {
      std::string arg = std::string(argv[++opt]);

      if (arg == "Generated") {
        tgttune = TargetTuneType::kGenerated;
      } else if (arg == "Tuned") {
        tgttune = TargetTuneType::kTuned;
      } else if (arg == "Tuned_BANFFPre") {
        tgttune = TargetTuneType::kTuned_BANFFPre;
      } else if (arg == "Tuned_BANFFPost") {
        tgttune = TargetTuneType::kTuned_BANFFPost;
      } else if (arg == "Mnv1Pi") {
        tgttune = TargetTuneType::kMnv1Pi;
      } else if (arg == "NonQE") {
        tgttune = TargetTuneType::kNonQE;
      } else {
        std::cout << "[ERROR]: Invalid option passed to --To " << arg
                  << std::endl;
        abort();
      }

    } else if (std::string(argv[opt]) == "-M") {
      bymode = true;
    }
    opt++;
  }
}

int main(int argc, char const *argv[]) {
  DeclareColors();

  handleOpts(argc, argv);

  gStyle->SetOptStat(false);
  ValidPlots(input_file, JustDoT2K);
}