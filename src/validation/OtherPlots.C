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

std::vector<TH1 *> Split(TH1 *in, bool splity = true) {
  if (in->GetDimension() != 2) {
    std::cout << "[ERROR]: Attempting to split a TH" << in->GetDimension()
              << ", can only apply to TH2s." << std::endl;
    abort();
  }

  TH2 *in2 = static_cast<TH2 *>(in);
  std::vector<TH1 *> rtnvect;
  if (splity) {
    for (int i = 0; i < in2->GetYaxis()->GetNbins(); ++i) {
      std::stringstream ss("");
      ss << in2->GetName() << "_px" << (i + 1);
      rtnvect.push_back(in2->ProjectionX(ss.str().c_str(), i + 1, i + 1, "e"));
      ss.str("");

      ss << in2->GetYaxis()->GetBinLowEdge(i + 1) << " < "
         << in2->GetYaxis()->GetTitle() << " < "
         << in2->GetYaxis()->GetBinUpEdge(i + 1);

      rtnvect.back()->SetTitle(ss.str().c_str());
    }
  } else {
    for (int i = 0; i < in2->GetXaxis()->GetNbins(); ++i) {
      std::stringstream ss("");
      ss << in2->GetName() << "_py" << (i + 1);
      rtnvect.push_back(in2->ProjectionY(ss.str().c_str(), i + 1, i + 1, "e"));
      ss.str("");

      ss << in2->GetXaxis()->GetBinLowEdge(i + 1) << " < "
         << in2->GetXaxis()->GetTitle() << " < "
         << in2->GetXaxis()->GetBinUpEdge(i + 1);

      rtnvect.back()->SetTitle(ss.str().c_str());
    }
  }
  return rtnvect;
}

double GetMaximumBinPlusError(TH1 *h) {

  double max = -std::numeric_limits<double>::max();

  for (int i = 0; i < h->GetXaxis()->GetNbins(); ++i) {
    max = std::max(max, h->GetBinContent(i + 1) + h->GetBinError(i + 1));
  }

  return max;
}

double GetMaximumBin(TH1 *h) {

  double max = -std::numeric_limits<double>::max();

  for (int i = 0; i < h->GetXaxis()->GetNbins(); ++i) {
    max = std::max(max, h->GetBinContent(i + 1));
  }

  return max;
}

struct hblob {
  std::unique_ptr<TH1> From;
  std::unique_ptr<TH1> Target;
  std::vector<std::unique_ptr<TH1>> ReWeights;

  bool is2D;
  Detector _det;

  std::string varname;
  selection sel;

  void Load(std::unique_ptr<TFile> &fin, Detector det, std::string varname,
            nuspecies nuspec, std::string tgtstr, selection sel,
            TargetTuneType tgttune = TargetTuneType::kTuned, int mode = 0) {
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

    std::string suffix = std::to_string(det) + "_" + varname + "_" +
                         std::to_string(nuspec) + "_" + tgtstr + "_" +
                         std::to_string(static_cast<int>(tgttune)) + "_" +
                         SelectionList[sel];

    if (From) {
      From->SetName(("FROM_" + suffix).c_str());

      std::cout << "read " << fromhist_name << " as " << From->GetName()
                << std::endl;
    }
    if (Target) {
      Target->SetName(("TARGET_" + suffix).c_str());

      std::cout << "read " << targethist_name << " as " << Target->GetName()
                << std::endl;
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
          (std::string("ReWeight_") + suffix + std::to_string(i++)).c_str());
      std::cout << "read " << n << " as " << ReWeights.back()->GetName()
                << std::endl;
    }
  }
};

void ValidPlots(std::string const &finname) {
  std::unique_ptr<TFile> fin(new TFile(finname.c_str()));
  if (fin->IsZombie()) {
    std::cout << "Failed to read " << finname << std::endl;
    return;
  }

  if (true) {
    for (std::string tarstr : {"CH", "H2O"}) {
      for (std::string varstr : {"q0", "Enu"}) {
        for (auto sel : {kCC0pi, kCC0pi_QE, kCC0pi_2p2h, kCCInc}) {

          std::string selstr = SelectionList[sel];

          hblob CH_Enu_numu;
          CH_Enu_numu.Load(fin, kT2K, varstr, kNuMu, tarstr, sel);
          hblob CH_Enu_numub;
          CH_Enu_numub.Load(fin, kT2K, varstr, kNuMub, tarstr, selection(sel));

          CH_Enu_numu.From->Divide(CH_Enu_numub.From.get());
          CH_Enu_numu.Target->Divide(CH_Enu_numub.Target.get());

          double max_gen = GetMaximumTH1s(
              std::vector<std::reference_wrapper<std::unique_ptr<TH1>>>{
                  CH_Enu_numu.From, CH_Enu_numu.Target});

          TCanvas *c1 = MakeCanvasTopLegend();
          TPad *p1 = MakeRatioTopPadTopLegend();
          p1->AppendPad();
          TPad *p2 = MakeRatioBottomPadTopLegend();
          p2->AppendPad();

          p1->cd();

          CH_Enu_numu.From->SetLineWidth(3);
          CH_Enu_numu.From->SetLineColor(SORNBrightWheel[1]);

          CH_Enu_numu.Target->SetLineWidth(3);
          CH_Enu_numu.Target->SetLineColor(SORNBrightWheel[2]);

          if (varstr == "Enu") {
            CH_Enu_numu.From->GetXaxis()->SetRangeUser(0, 2);
          }

          CH_Enu_numu.From->GetYaxis()->SetRangeUser(0, max_gen * 1.01);
          CH_Enu_numu.From->GetYaxis()->SetTitle(
              "#sigma_{#nu}/#sigma_{#bar{#nu}}");
          CH_Enu_numu.From->SetTitle("");

          CH_Enu_numu.From->GetYaxis()->SetNdivisions(505);
          CH_Enu_numu.From->GetYaxis()->SetLabelSize(0.07);
          CH_Enu_numu.From->GetYaxis()->SetLabelOffset(0.01);
          CH_Enu_numu.From->GetYaxis()->SetTitleSize(0.07);
          CH_Enu_numu.From->GetYaxis()->SetTitleOffset(1);
          HideAxis(CH_Enu_numu.From->GetXaxis());

          CH_Enu_numu.From->DrawClone("HIST");
          CH_Enu_numu.Target->DrawClone("HISTSAME");

          c1->cd();
          TLegend *leg = new TLegend(0.125, 0.81, 0.925, 1);
          leg->SetNColumns(2);
          leg->SetTextSize(0.06);
          leg->SetBorderSize(0);
          leg->SetFillStyle(4001);

          leg->AddEntry(CH_Enu_numu.From.get(), "T2K Post. ND", "l");
          leg->AddEntry(CH_Enu_numu.Target.get(), "NOvA2020", "l");
          leg->Draw();

          p2->cd();

          CH_Enu_numu.From->Divide(CH_Enu_numu.Target.get());
          CH_Enu_numu.From->GetYaxis()->SetRangeUser(0.5, 2);

          CH_Enu_numu.From->GetYaxis()->SetNdivisions(505);
          CH_Enu_numu.From->GetYaxis()->SetLabelSize(0.07 * 2);
          CH_Enu_numu.From->GetYaxis()->SetLabelOffset(0.01);
          CH_Enu_numu.From->GetYaxis()->SetTitleSize(0.07 * 2);
          CH_Enu_numu.From->GetYaxis()->SetTitleOffset(1. / 2.);

          CH_Enu_numu.From->GetXaxis()->SetNdivisions(505);
          CH_Enu_numu.From->GetXaxis()->SetLabelSize(0.07 * 2);
          CH_Enu_numu.From->GetXaxis()->SetLabelOffset(0.01);
          CH_Enu_numu.From->GetXaxis()->SetTitleSize(0.07 * 2);
          CH_Enu_numu.From->GetXaxis()->SetTitleOffset(1.1);

          CH_Enu_numu.From->GetYaxis()->SetTitle("NEUT/GENIE");
          CH_Enu_numu.From->SetTitle("");
          CH_Enu_numu.From->DrawClone("HIST");

          c1->cd();

          TLatex ttl;
          ttl.SetTextSize(0.05);
          ttl.SetTextAlign(32);
          ttl.DrawLatexNDC(0.95, 0.775, selstr.c_str());

          c1->Print((std::string("double_ratio_nunubar_") + tarstr + "_" +
                     selstr + "_" + varstr + ".pdf")
                        .c_str());
        }
      }
    }

    for (std::string varstr : {"q0", "Enu"}) {
      for (auto sel : {kCC0pi, kCC0pi_QE, kCC0pi_2p2h, kCCInc}) {

        std::string selstr = SelectionList[sel];

        hblob CH_Enu_numu;
        CH_Enu_numu.Load(fin, kT2K, varstr, kNuMu, "H2O", selection(sel));
        hblob CH_Enu_numub;
        CH_Enu_numub.Load(fin, kT2K, varstr, kNuMu, "CH", selection(sel));

        CH_Enu_numu.From->Divide(CH_Enu_numub.From.get());
        CH_Enu_numu.Target->Divide(CH_Enu_numub.Target.get());

        double max_gen = GetMaximumTH1s(
            std::vector<std::reference_wrapper<std::unique_ptr<TH1>>>{
                CH_Enu_numu.From, CH_Enu_numu.Target});

        TCanvas *c1 = MakeCanvasTopLegend();
        TPad *p1 = MakeRatioTopPadTopLegend();
        p1->AppendPad();
        TPad *p2 = MakeRatioBottomPadTopLegend();
        p2->AppendPad();

        p1->cd();

        CH_Enu_numu.From->SetLineWidth(3);
        CH_Enu_numu.From->SetLineColor(SORNBrightWheel[1]);

        CH_Enu_numu.Target->SetLineWidth(3);
        CH_Enu_numu.Target->SetLineColor(SORNBrightWheel[2]);

        CH_Enu_numu.From->GetYaxis()->SetRangeUser(0, max_gen * 1.01);
        CH_Enu_numu.From->GetYaxis()->SetTitle("#sigma_{H2O}/#sigma_{CH}");
        CH_Enu_numu.From->SetTitle("");

        CH_Enu_numu.From->GetYaxis()->SetNdivisions(505);
        CH_Enu_numu.From->GetYaxis()->SetLabelSize(0.07);
        CH_Enu_numu.From->GetYaxis()->SetLabelOffset(0.01);
        CH_Enu_numu.From->GetYaxis()->SetTitleSize(0.07);
        CH_Enu_numu.From->GetYaxis()->SetTitleOffset(1);
        HideAxis(CH_Enu_numu.From->GetXaxis());

        CH_Enu_numu.From->DrawClone("HIST");
        CH_Enu_numu.Target->DrawClone("HISTSAME");

        c1->cd();
        TLegend *leg = new TLegend(0.125, 0.81, 0.925, 1);
        leg->SetNColumns(2);
        leg->SetTextSize(0.06);
        leg->SetBorderSize(0);
        leg->SetFillStyle(4001);

        leg->AddEntry(CH_Enu_numu.From.get(), "T2K Post. ND", "l");
        leg->AddEntry(CH_Enu_numu.Target.get(), "NOvA2020", "l");
        leg->Draw();

        p2->cd();

        CH_Enu_numu.From->Divide(CH_Enu_numu.Target.get());
        CH_Enu_numu.From->GetYaxis()->SetRangeUser(0.5, 2);

        CH_Enu_numu.From->GetYaxis()->SetNdivisions(505);
        CH_Enu_numu.From->GetYaxis()->SetLabelSize(0.07 * 2);
        CH_Enu_numu.From->GetYaxis()->SetLabelOffset(0.01);
        CH_Enu_numu.From->GetYaxis()->SetTitleSize(0.07 * 2);
        CH_Enu_numu.From->GetYaxis()->SetTitleOffset(1. / 2.);

        CH_Enu_numu.From->GetXaxis()->SetNdivisions(505);
        CH_Enu_numu.From->GetXaxis()->SetLabelSize(0.07 * 2);
        CH_Enu_numu.From->GetXaxis()->SetLabelOffset(0.01);
        CH_Enu_numu.From->GetXaxis()->SetTitleSize(0.07 * 2);
        CH_Enu_numu.From->GetXaxis()->SetTitleOffset(1.1);

        CH_Enu_numu.From->GetYaxis()->SetTitle("NEUT/GENIE");
        CH_Enu_numu.From->SetTitle("");
        CH_Enu_numu.From->DrawClone("HIST");

        c1->cd();

        TLatex ttl;
        ttl.SetTextSize(0.05);
        ttl.SetTextAlign(32);
        ttl.DrawLatexNDC(0.95, 0.775, selstr.c_str());

        c1->Print((std::string("double_ratio_CHH2O_") + selstr + "_" + varstr +
                   ".pdf")
                      .c_str());
      }
    }
  }

  int colorwheel[6] = {SORNVibrantWheel[1], SORNVibrantWheel[3], kBlack,
                       SORNVibrantWheel[2]};

  { // How many panes do we need
    int nslices = 8;
    int first_slice = 1;

    // don't get more than 3 wide
    int nx = 3;
    int ny = 3;

    std::cout << "[INFO]: Build " << nx << "x" << ny << " pads for " << nslices
              << std::endl;

    if ((nx * ny) < nslices) {
      std::cout << "[ERROR]: With " << nslices
                << " total, we calculated that we need nx = " << nx
                << ", ny = " << ny << ", total = " << (nx * ny) << " panes."
                << std::endl;
    }

    double ybuffer = 0;

    if ((nx * ny) ==
        nslices) { // if we have a full house, put the legend at the top
      ybuffer = 0.15;
    }

    double padxwidth = 1.0 / double(nx);
    double padywidth = (1.0 - ybuffer) / double(ny);

    TCanvas c1("c1", "", nx * 400, (ny * 400) * (1 + ybuffer));

    c1.Print("T2K_ERecQEBiasSlices.pdf[");

    for (auto sel : {kCCInc, kCC0pi, kCC0pi_QE, kCC0pi_2p2h, kCC1cpi, kCC1pi0,
                     kCCmultipi, kCCOther}) {

      c1.Clear();

      hblob CH_numu_EnuERecQEBias;
      CH_numu_EnuERecQEBias.Load(fin, kT2K, "EnuERecQEBias", kNuMu, "CH", sel);

      hblob CH_numu_EnuERecQEBias_Mnv1Pi;
      CH_numu_EnuERecQEBias_Mnv1Pi.Load(fin, kT2K, "EnuERecQEBias", kNuMu, "CH",
                                        sel, TargetTuneType::kMnv1Pi);

      hblob CH_numu_EnuERecQEBias_NonQE;
      CH_numu_EnuERecQEBias_NonQE.Load(fin, kT2K, "EnuERecQEBias", kNuMu,
      "CH",
                                       sel, TargetTuneType::kNonQE);

      auto splits_BANFF = Split(CH_numu_EnuERecQEBias.From.get(), false);
      auto splits_Mnv1Pi =
          Split(CH_numu_EnuERecQEBias_Mnv1Pi.Target.get(), false);
      auto splits_NonQE =
          Split(CH_numu_EnuERecQEBias_NonQE.Target.get(), false);
      auto splits_NOvA = Split(CH_numu_EnuERecQEBias.Target.get(), false);

      std::stringstream ss("");
      for (int i = first_slice; i < (nslices + first_slice); ++i) {
        int ix = (i - first_slice) % nx;
        int iy = (i - first_slice) / nx;

        double px1 = ix * padxwidth;
        double px2 = (ix + 1) * padxwidth;
        double py1 = (1.0 - ybuffer) - ((iy + 1) * padywidth);
        double py2 = (1.0 - ybuffer) - (iy * padywidth);

        ss.str("");
        ss << "p" << i;
        TPad *p = new TPad(ss.str().c_str(), "", px1, py1, px2, py2);
        p->SetBottomMargin(0.15);
        p->SetLeftMargin(0.15);
        p->AppendPad();
        p->cd();

        splits_NOvA[i]->Divide(splits_BANFF[i]);
        splits_Mnv1Pi[i]->Divide(splits_BANFF[i]);
        splits_NonQE[i]->Divide(splits_BANFF[i]);
        splits_BANFF[i]->Divide(splits_BANFF[i]);

        splits_BANFF[i]->GetYaxis()->SetRangeUser(0.5, 2);

        splits_BANFF[i]->GetXaxis()->SetNdivisions(505);
        splits_BANFF[i]->GetYaxis()->SetNdivisions(505);

        splits_BANFF[i]->GetYaxis()->SetTitle("Relative rate to T2K Post. ND");
        splits_BANFF[i]->GetYaxis()->SetTitleSize(0.06);
        splits_BANFF[i]->GetYaxis()->SetTitleOffset(1.2);
        splits_BANFF[i]->GetYaxis()->SetLabelSize(0.06);

        splits_BANFF[i]->GetXaxis()->SetTitleSize(0.06);
        splits_BANFF[i]->GetXaxis()->SetTitleOffset(1);
        splits_BANFF[i]->GetXaxis()->SetLabelSize(0.06);

        splits_BANFF[i]->SetLineColor(colorwheel[2]);
        splits_BANFF[i]->SetLineWidth(2);
        splits_BANFF[i]->SetLineStyle(1);
        splits_BANFF[i]->SetMarkerSize(0);
        splits_BANFF[i]->SetMarkerStyle(0);
        splits_BANFF[i]->SetMarkerColorAlpha(colorwheel[2], 0);

        splits_NOvA[i]->SetLineColor(colorwheel[0]);
        splits_NOvA[i]->SetLineWidth(2);
        splits_NOvA[i]->SetLineStyle(1);
        splits_NOvA[i]->SetMarkerSize(0);
        splits_NOvA[i]->SetMarkerStyle(0);
        splits_NOvA[i]->SetMarkerColorAlpha(colorwheel[0], 0);

        splits_Mnv1Pi[i]->SetLineColor(colorwheel[1]);
        splits_Mnv1Pi[i]->SetLineWidth(2);
        splits_Mnv1Pi[i]->SetLineStyle(1);
        splits_Mnv1Pi[i]->SetMarkerSize(0);
        splits_Mnv1Pi[i]->SetMarkerStyle(0);
        splits_Mnv1Pi[i]->SetMarkerColorAlpha(colorwheel[1], 0);

        splits_NonQE[i]->SetLineColor(colorwheel[3]);
        splits_NonQE[i]->SetLineWidth(2);
        splits_NonQE[i]->SetLineStyle(1);
        splits_NonQE[i]->SetMarkerSize(0);
        splits_NonQE[i]->SetMarkerStyle(0);
        splits_NonQE[i]->SetMarkerColorAlpha(colorwheel[3], 0);

        splits_BANFF[i]->DrawClone("HIST");
        splits_NOvA[i]->DrawClone("HISTSAME");
        splits_Mnv1Pi[i]->DrawClone("HISTSAME");
        splits_NonQE[i]->DrawClone("HISTSAME");

        c1.cd();
      }

      TLegend *leg = nullptr;

      if (ybuffer > 0) {
        leg = new TLegend(0.1, 0.85, 0.9, 1);
        leg->SetTextSize(0.03);
      } else { // Put the legend in the place of the last pane(s)
        leg = new TLegend((nslices % nx) * padxwidth, 1 - (ny * padywidth),
                          nx * padxwidth, 1 - ((ny - 1) * padywidth * 1.1));
        leg->SetTextSize(0.03 * float((nx * ny) - nslices));
      }

      leg->SetBorderSize(0);
      leg->SetFillStyle(0);
      leg->SetNColumns(1);
      leg->AddEntry(splits_BANFF[first_slice], "T2K Post. ND", "lp");
      leg->AddEntry(splits_NOvA[first_slice], "NOvA NDTune", "lp");
      leg->AddEntry(splits_Mnv1Pi[first_slice], "T2K Pre. + Mnv1Pi", "lp");
      leg->AddEntry(splits_NonQE[first_slice], "T2K NonQE", "lp");

      leg->Draw();

      std::string selstr = SelectionList[sel];

      TLatex ttl;
      ttl.SetTextSize(0.05);
      ttl.SetTextAlign(32);
      ttl.DrawLatexNDC(0.95, 0.3, selstr.c_str());

      c1.Print("T2K_ERecQEBiasSlices.pdf");
    }
    c1.Print("T2K_ERecQEBiasSlices.pdf]");
  }

  { // How many panes do we need
    int nslices = 8;
    int first_slice = 1;

    // don't get more than 3 wide
    int nx = 3;
    int ny = 3;

    std::cout << "[INFO]: Build " << nx << "x" << ny << " pads for " << nslices
              << std::endl;

    if ((nx * ny) < nslices) {
      std::cout << "[ERROR]: With " << nslices
                << " total, we calculated that we need nx = " << nx
                << ", ny = " << ny << ", total = " << (nx * ny) << " panes."
                << std::endl;
    }

    double ybuffer = 0;

    if ((nx * ny) ==
        nslices) { // if we have a full house, put the legend at the top
      ybuffer = 0.15;
    }

    double padxwidth = 1.0 / double(nx);
    double padywidth = (1.0 - ybuffer) / double(ny);

    TCanvas c1("c1", "", nx * 400, (ny * 400) * (1 + ybuffer));

    c1.Print("NOvA_ERecAvBiasSlices.pdf[");

    for (auto sel : {kCCInc, kCC0pi, kCC0pi_QE, kCC0pi_2p2h, kCC1cpi, kCC1pi0,
                     kCCmultipi, kCCOther}) {

      c1.Clear();

      hblob CH_numu_EnuERecAvBias;
      CH_numu_EnuERecAvBias.Load(fin, kNOvAND, "EnuERecAvBias", kNuMu, "CH",
                                 selection(sel));

      hblob CH_numu_EnuERecAvBias_Mnv1Pi;
      CH_numu_EnuERecAvBias_Mnv1Pi.Load(fin, kNOvAND, "EnuERecAvBias", kNuMu,
                                        "CH", sel, TargetTuneType::kMnv1Pi);

      hblob CH_numu_EnuERecAvBias_NonQE;
      CH_numu_EnuERecAvBias_NonQE.Load(fin, kNOvAND, "EnuERecAvBias", kNuMu,
                                       "CH", sel, TargetTuneType::kNonQE);

      auto splits_BANFF = Split(CH_numu_EnuERecAvBias.Target.get(), false);
      auto splits_Mnv1Pi =
          Split(CH_numu_EnuERecAvBias_Mnv1Pi.Target.get(), false);
      auto splits_NonQE =
          Split(CH_numu_EnuERecAvBias_NonQE.Target.get(), false);
      auto splits_NOvA = Split(CH_numu_EnuERecAvBias.From.get(), false);

      std::stringstream ss("");
      for (int i = first_slice; i < (nslices + first_slice); ++i) {
        int ix = (i - first_slice) % nx;
        int iy = (i - first_slice) / nx;

        double px1 = ix * padxwidth;
        double px2 = (ix + 1) * padxwidth;
        double py1 = (1.0 - ybuffer) - ((iy + 1) * padywidth);
        double py2 = (1.0 - ybuffer) - (iy * padywidth);

        ss.str("");
        ss << "p" << i;
        TPad *p = new TPad(ss.str().c_str(), "", px1, py1, px2, py2);
        p->SetBottomMargin(0.15);
        p->SetLeftMargin(0.15);
        p->AppendPad();
        p->cd();

        splits_Mnv1Pi[i]->Divide(splits_NOvA[i]);
        splits_BANFF[i]->Divide(splits_NOvA[i]);
        splits_NonQE[i]->Divide(splits_NOvA[i]);
        splits_NOvA[i]->Divide(splits_NOvA[i]);

        splits_NOvA[i]->GetYaxis()->SetRangeUser(0.5, 1.5);
        splits_NOvA[i]->GetXaxis()->SetRangeUser(-0.3, 0.05);

        splits_NOvA[i]->GetXaxis()->SetNdivisions(505);
        splits_NOvA[i]->GetYaxis()->SetNdivisions(505);

        splits_NOvA[i]->GetYaxis()->SetTitle("Relative rate to NOvA NDTune");
        splits_NOvA[i]->GetYaxis()->SetTitleSize(0.06);
        splits_NOvA[i]->GetYaxis()->SetTitleOffset(1.2);
        splits_NOvA[i]->GetYaxis()->SetLabelSize(0.06);

        splits_NOvA[i]->GetXaxis()->SetTitleSize(0.06);
        splits_NOvA[i]->GetXaxis()->SetTitleOffset(1);
        splits_NOvA[i]->GetXaxis()->SetLabelSize(0.06);

        splits_NOvA[i]->SetLineColor(colorwheel[0]);
        splits_NOvA[i]->SetLineWidth(2);
        splits_NOvA[i]->SetLineStyle(1);
        splits_NOvA[i]->SetMarkerSize(0);
        splits_NOvA[i]->SetMarkerStyle(0);
        splits_NOvA[i]->SetMarkerColorAlpha(colorwheel[0], 0);

        splits_BANFF[i]->SetLineColor(colorwheel[2]);
        splits_BANFF[i]->SetLineWidth(2);
        splits_BANFF[i]->SetLineStyle(1);
        splits_BANFF[i]->SetMarkerSize(0);
        splits_BANFF[i]->SetMarkerStyle(0);
        splits_BANFF[i]->SetMarkerColorAlpha(colorwheel[2], 0);

        splits_Mnv1Pi[i]->SetLineColor(colorwheel[1]);
        splits_Mnv1Pi[i]->SetLineWidth(2);
        splits_Mnv1Pi[i]->SetLineStyle(1);
        splits_Mnv1Pi[i]->SetMarkerSize(0);
        splits_Mnv1Pi[i]->SetMarkerStyle(0);
        splits_Mnv1Pi[i]->SetMarkerColorAlpha(colorwheel[1], 0);

        splits_NonQE[i]->SetLineColor(colorwheel[3]);
        splits_NonQE[i]->SetLineWidth(2);
        splits_NonQE[i]->SetLineStyle(1);
        splits_NonQE[i]->SetMarkerSize(0);
        splits_NonQE[i]->SetMarkerStyle(0);
        splits_NonQE[i]->SetMarkerColorAlpha(colorwheel[3], 0);

        if (splits_BANFF[i]->Integral() > 0) {

          splits_NOvA[i]->DrawClone("HIST");
          splits_BANFF[i]->DrawClone("HISTSAME");
          splits_Mnv1Pi[i]->DrawClone("HISTSAME");
          splits_NonQE[i]->DrawClone("HISTSAME");
        }

        c1.cd();
      }

      TLegend *leg = nullptr;

      if (ybuffer > 0) {
        leg = new TLegend(0.1, 0.85, 0.9, 1);
        leg->SetTextSize(0.03);
      } else { // Put the legend in the place of the last pane(s)
        leg = new TLegend((nslices % nx) * padxwidth, 1 - (ny * padywidth),
                          nx * padxwidth, 1 - ((ny - 1) * padywidth * 1.1));
        leg->SetTextSize(0.03 * float((nx * ny) - nslices));
      }

      leg->SetBorderSize(0);
      leg->SetFillStyle(0);
      leg->SetNColumns(1);
      leg->AddEntry(splits_BANFF[first_slice], "T2K Post. ND", "lp");
      leg->AddEntry(splits_NOvA[first_slice], "NOvA NDTune", "lp");
      leg->AddEntry(splits_Mnv1Pi[first_slice], "T2K Pre. + Mnv1Pi", "lp");
      leg->AddEntry(splits_NonQE[first_slice], "T2K NonQE", "lp");

      leg->Draw();

      std::string selstr = SelectionList[sel];

      TLatex ttl;
      ttl.SetTextSize(0.05);
      ttl.SetTextAlign(32);
      ttl.DrawLatexNDC(0.95, 0.3, selstr.c_str());

      c1.Print("NOvA_ERecAvBiasSlices.pdf");
    }
    c1.Print("NOvA_ERecAvBiasSlices.pdf]");
  }
}

int main(int argc, char const *argv[]) {
  DeclareColors();
  gStyle->SetOptStat(false);
  ValidPlots(argv[1]);
}