#include "T2KNOvA/FakeDataHelper.hxx"
#include "TLatex.h"
#include "colordef.h"
#include "plotutils.h"

#include "TExec.h"

using namespace t2knova;

bool const kT2K = true;
bool const kNOvAND = false;
bool bymode = false;
bool doosc = false;

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

  void Load(std::unique_ptr<TFile> &fin, bool t2kbase, bool istuned,
            std::string varname, nuspecies nuspec, std::string tgtstr,
            selection sel, int mode = 0) {

    std::string mode_str = "";
    if (mode != 0) {
      mode_str =
          (mode < 0 ? "_Mode_m" : "_Mode_") + std::to_string(std::abs(mode));
    }

    std::string sel_str = "";
    if (sel != t2knova::kNoPrimarySel) {
      sel_str = std::string("_") + SelectionList[sel];
    }

    std::string frompath =
        (t2kbase ? std::string("ND280/") + (istuned ? "T2KNDTune" : "NEUT")
                 : std::string("NOvAND/") + (istuned ? "NOvATune" : "GENIE")) +
        +"/" + tgtstr + "/" + all_nuspecies[nuspec] + "/" + varname + sel_str +
        mode_str;

    std::string targetpath =
        (t2kbase ? std::string("ND280/") + (istuned ? "NOvATune" : "GENIE")
                 : std::string("NOvAND/") + (istuned ? "T2KNDTune" : "NEUT")) +
        +"/" + tgtstr + "/" + all_nuspecies[nuspec] + "/" + varname + sel_str +
        mode_str;

    From = GetTH1(fin, frompath, false);

    Target = GetTH1(fin, targetpath, false);

    if (From) {
      From->SetName("FROM");
    } else {
      std::cout << "[ERROR]: " << frompath << std::endl;
      abort();
    }
    if (Target) {
      Target->SetName("TARGET");
    } else {
      std::cout << "[ERROR]: " << targetpath << std::endl;
      abort();
    }
  }
};

void ValidPlots(std::string const &finname) {
  std::unique_ptr<TFile> fin(new TFile(finname.c_str()));
  if (fin->IsZombie()) {
    std::cout << "Failed to read " << finname << std::endl;
    return;
  }

  if (false) {
    for (std::string tarstr : {"CH", "H2O"}) {
      for (std::string varstr : {"q0", "Enu"}) {
        for (auto sel : {kCC0pi, kCC0pi_QE, kCC0pi_2p2h, kCCInc}) {

          std::string selstr = SelectionList[sel];

          hblob CH_Enu_numu;
          CH_Enu_numu.Load(fin, kT2K, true, varstr, kNuMu, tarstr,
                           selection(sel));
          hblob CH_Enu_numub;
          CH_Enu_numub.Load(fin, kT2K, true, varstr, kNuMub, tarstr,
                            selection(sel));
          hblob CH_Enu_numu_untuned;
          CH_Enu_numu_untuned.Load(fin, kT2K, false, varstr, kNuMu, tarstr,
                                   selection(sel));
          hblob CH_Enu_numub_untuned;
          CH_Enu_numub_untuned.Load(fin, kT2K, false, varstr, kNuMub, tarstr,
                                    selection(sel));

          CH_Enu_numu.From->Divide(CH_Enu_numub.From.get());
          CH_Enu_numu.Target->Divide(CH_Enu_numub.Target.get());
          CH_Enu_numu_untuned.From->Divide(CH_Enu_numub_untuned.From.get());
          CH_Enu_numu_untuned.Target->Divide(CH_Enu_numub_untuned.Target.get());

          double max_gen = GetMaximumTH1s(
              std::vector<std::reference_wrapper<std::unique_ptr<TH1>>>{
                  CH_Enu_numu.From, CH_Enu_numu.Target,
                  CH_Enu_numu_untuned.From, CH_Enu_numu_untuned.Target});

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

          CH_Enu_numu_untuned.From->SetLineStyle(2);
          CH_Enu_numu_untuned.From->SetLineWidth(3);
          CH_Enu_numu_untuned.From->SetLineColor(SORNBrightWheel[3]);

          CH_Enu_numu_untuned.Target->SetLineStyle(2);
          CH_Enu_numu_untuned.Target->SetLineWidth(3);
          CH_Enu_numu_untuned.Target->SetLineColor(SORNBrightWheel[4]);

          if (varstr == "Enu") {
            CH_Enu_numu.From->GetXaxis()->SetRangeUser(0, 2);
          }

          CH_Enu_numu.From->GetYaxis()->SetRangeUser(0, max_gen * 1.01);
          CH_Enu_numu.From->GetYaxis()->SetTitle(
              "#sigma_{#nu}/#sigma_{#bar{#nu}}");
          // CH_Enu_numu.From->GetYaxis()->SetTitle("#sigma_{#nu}/#sigma_{#bar{#nu}}");
          CH_Enu_numu.From->SetTitle("");

          CH_Enu_numu.From->GetYaxis()->SetNdivisions(505);
          CH_Enu_numu.From->GetYaxis()->SetLabelSize(0.07);
          CH_Enu_numu.From->GetYaxis()->SetLabelOffset(0.01);
          CH_Enu_numu.From->GetYaxis()->SetTitleSize(0.07);
          CH_Enu_numu.From->GetYaxis()->SetTitleOffset(1);
          HideAxis(CH_Enu_numu.From->GetXaxis());

          CH_Enu_numu.From->DrawClone("HIST");
          CH_Enu_numu.Target->DrawClone("HISTSAME");
          CH_Enu_numu_untuned.From->DrawClone("HISTSAME");
          CH_Enu_numu_untuned.Target->DrawClone("HISTSAME");

          c1->cd();
          TLegend *leg = new TLegend(0.125, 0.81, 0.925, 1);
          leg->SetNColumns(2);
          leg->SetTextSize(0.06);
          leg->SetBorderSize(0);
          leg->SetFillStyle(4001);

          leg->AddEntry(CH_Enu_numu.From.get(), "BANFF", "l");
          leg->AddEntry(CH_Enu_numu.Target.get(), "NOvA2020", "l");
          leg->AddEntry(CH_Enu_numu_untuned.From.get(), "NEUT", "l");
          leg->AddEntry(CH_Enu_numu_untuned.Target.get(), "GENIE", "l");
          leg->Draw();

          p2->cd();

          CH_Enu_numu.From->Divide(CH_Enu_numu.Target.get());
          CH_Enu_numu_untuned.From->Divide(CH_Enu_numu_untuned.Target.get());
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
          CH_Enu_numu_untuned.From->DrawClone("HISTSAME");

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
        CH_Enu_numu.Load(fin, true, true, varstr, kNuMu, "H2O", selection(sel));
        hblob CH_Enu_numub;
        CH_Enu_numub.Load(fin, true, true, varstr, kNuMu, "CH", selection(sel));
        hblob CH_Enu_numu_untuned;
        CH_Enu_numu_untuned.Load(fin, true, false, varstr, kNuMu, "H2O",
                                 selection(sel));
        hblob CH_Enu_numub_untuned;
        CH_Enu_numub_untuned.Load(fin, true, false, varstr, kNuMu, "CH",
                                  selection(sel));

        CH_Enu_numu.From->Divide(CH_Enu_numub.From.get());
        CH_Enu_numu.Target->Divide(CH_Enu_numub.Target.get());
        CH_Enu_numu_untuned.From->Divide(CH_Enu_numub_untuned.From.get());
        CH_Enu_numu_untuned.Target->Divide(CH_Enu_numub_untuned.Target.get());

        double max_gen = GetMaximumTH1s(
            std::vector<std::reference_wrapper<std::unique_ptr<TH1>>>{
                CH_Enu_numu.From, CH_Enu_numu.Target, CH_Enu_numu_untuned.From,
                CH_Enu_numu_untuned.Target});

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

        CH_Enu_numu_untuned.From->SetLineStyle(2);
        CH_Enu_numu_untuned.From->SetLineWidth(3);
        CH_Enu_numu_untuned.From->SetLineColor(SORNBrightWheel[3]);

        CH_Enu_numu_untuned.Target->SetLineStyle(2);
        CH_Enu_numu_untuned.Target->SetLineWidth(3);
        CH_Enu_numu_untuned.Target->SetLineColor(SORNBrightWheel[4]);

        CH_Enu_numu.From->GetYaxis()->SetRangeUser(0, max_gen * 1.01);
        // CH_Enu_numu.From->GetYaxis()->SetTitle(
        //     "#sigma_{#nu}/#sigma_{#bar{#nu}}");
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
        CH_Enu_numu_untuned.From->DrawClone("HISTSAME");
        CH_Enu_numu_untuned.Target->DrawClone("HISTSAME");

        c1->cd();
        TLegend *leg = new TLegend(0.125, 0.81, 0.925, 1);
        leg->SetNColumns(2);
        leg->SetTextSize(0.06);
        leg->SetBorderSize(0);
        leg->SetFillStyle(4001);

        leg->AddEntry(CH_Enu_numu.From.get(), "BANFF", "l");
        leg->AddEntry(CH_Enu_numu.Target.get(), "NOvA2020", "l");
        leg->AddEntry(CH_Enu_numu_untuned.From.get(), "NEUT", "l");
        leg->AddEntry(CH_Enu_numu_untuned.Target.get(), "GENIE", "l");
        leg->Draw();

        p2->cd();

        CH_Enu_numu.From->Divide(CH_Enu_numu.Target.get());
        CH_Enu_numu_untuned.From->Divide(CH_Enu_numu_untuned.Target.get());
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
        CH_Enu_numu_untuned.From->DrawClone("HISTSAME");

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

  int colorwheel[6] = {SORNVibrantWheel[1],SORNVibrantWheel[3],kBlack,SORNVibrantWheel[2]};

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

    for (int sel : {kCCInc, kCC0pi, kCC0pi_QE, kCC0pi_2p2h, kCC1cpi, kCC1pi0,
                    kCCmultipi, kCCOther}) {

      c1.Clear();

      hblob CH_numu_EnuERecQEBias;
      CH_numu_EnuERecQEBias.Load(fin, kT2K, true, "EnuERecQEBias", kNuMu, "CH",
                                 selection(sel));

      hblob CH_numu_EnuERecQEBias_untuned;
      CH_numu_EnuERecQEBias_untuned.Load(fin, kT2K, false, "EnuERecQEBias",
                                         kNuMu, "CH", selection(sel));

      auto splits_BANFF = Split(CH_numu_EnuERecQEBias.From.get(), false);
      auto splits_NOvA = Split(CH_numu_EnuERecQEBias.Target.get(), false);

      CH_numu_EnuERecQEBias_untuned.From->SetName(
          "CH_numu_EnuERecQEBias_untuned_Target");
      CH_numu_EnuERecQEBias_untuned.Target->SetName(
          "CH_numu_EnuERecQEBias_untuned_From");

      auto splits_BANFF_untuned =
          Split(CH_numu_EnuERecQEBias_untuned.Target.get(), false);
      auto splits_NOvA_untuned =
          Split(CH_numu_EnuERecQEBias_untuned.From.get(), false);

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

        splits_BANFF[i]->Scale(1.0 / GetMaximumBin(splits_BANFF[i]));
        splits_NOvA[i]->Scale(1.0 / GetMaximumBin(splits_NOvA[i]));
        splits_BANFF_untuned[i]->Scale(1.0 /
                                       GetMaximumBin(splits_BANFF_untuned[i]));
        splits_NOvA_untuned[i]->Scale(1.0 /
                                      GetMaximumBin(splits_NOvA_untuned[i]));

        splits_BANFF[i]->GetYaxis()->SetRangeUser(0, 1.01);

        splits_BANFF[i]->GetXaxis()->SetNdivisions(505);
        splits_BANFF[i]->GetYaxis()->SetNdivisions(505);

        splits_BANFF[i]->GetYaxis()->SetTitle("Relative rate");
        splits_BANFF[i]->GetYaxis()->SetTitleSize(0.06);
        splits_BANFF[i]->GetYaxis()->SetTitleOffset(0.8);
        splits_BANFF[i]->GetYaxis()->SetLabelSize(0.06);

        splits_BANFF[i]->GetXaxis()->SetTitleSize(0.06);
        splits_BANFF[i]->GetXaxis()->SetTitleOffset(1);
        splits_BANFF[i]->GetXaxis()->SetLabelSize(0.06);

        splits_BANFF[i]->SetLineColor(colorwheel[0]);
        splits_BANFF[i]->SetLineWidth(2);
        splits_BANFF[i]->SetLineStyle(1);
        splits_BANFF[i]->SetMarkerSize(0);
        splits_BANFF[i]->SetMarkerStyle(0);
        splits_BANFF[i]->SetMarkerColorAlpha(colorwheel[0], 0);

        splits_BANFF_untuned[i]->SetLineColor(colorwheel[1]);
        splits_BANFF_untuned[i]->SetLineWidth(2);
        splits_BANFF_untuned[i]->SetLineStyle(2);
        splits_BANFF_untuned[i]->SetMarkerSize(0);
        splits_BANFF_untuned[i]->SetMarkerStyle(0);
        splits_BANFF_untuned[i]->SetMarkerColorAlpha(colorwheel[1], 0);

        splits_NOvA[i]->SetLineColor(colorwheel[2]);
        splits_NOvA[i]->SetLineWidth(2);
        splits_NOvA[i]->SetLineStyle(1);
        splits_NOvA[i]->SetMarkerSize(0);
        splits_NOvA[i]->SetMarkerStyle(0);
        splits_NOvA[i]->SetMarkerColorAlpha(colorwheel[2], 0);

        splits_NOvA_untuned[i]->SetLineColor(colorwheel[3]);
        splits_NOvA_untuned[i]->SetLineWidth(2);
        splits_NOvA_untuned[i]->SetLineStyle(2);
        splits_NOvA_untuned[i]->SetMarkerSize(0);
        splits_NOvA_untuned[i]->SetMarkerStyle(0);
        splits_NOvA_untuned[i]->SetMarkerColorAlpha(colorwheel[1], 0);

        splits_BANFF[i]->DrawClone("HIST");
        splits_NOvA[i]->DrawClone("HISTSAME");
        splits_BANFF_untuned[i]->DrawClone("HISTSAME");
        splits_NOvA_untuned[i]->DrawClone("HISTSAME");

        c1.cd();
      }

      TLegend *leg = nullptr;

      if (ybuffer > 0) {
        leg = new TLegend(0.1, 0.85, 0.9, 1);
        leg->SetTextSize(0.04);
      } else { // Put the legend in the place of the last pane(s)
        leg = new TLegend((nslices % nx) * padxwidth, 1 - (ny * padywidth),
                          nx * padxwidth, 1 - ((ny - 1) * padywidth));
        leg->SetTextSize(0.04 * float((nx * ny) - nslices));
      }

      leg->SetBorderSize(0);
      leg->SetFillStyle(0);
      leg->SetNColumns(2);
      leg->AddEntry(splits_BANFF[first_slice], "BANFF", "lp");
      leg->AddEntry(splits_BANFF_untuned[first_slice], "NEUT", "lp");
      leg->AddEntry(splits_NOvA[first_slice], "NOvA", "lp");
      leg->AddEntry(splits_NOvA_untuned[first_slice], "GENIE", "lp");

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

    for (int sel : {kCCInc, kCC0pi, kCC0pi_QE, kCC0pi_2p2h, kCC1cpi, kCC1pi0,
                    kCCmultipi, kCCOther}) {

      c1.Clear();

      hblob CH_numu_EnuERecAvBias;
      CH_numu_EnuERecAvBias.Load(fin, kNOvAND, true, "EnuERecAvBias", kNuMu,
                                 "CH", selection(sel));

      hblob CH_numu_EnuERecAvBias_untuned;
      CH_numu_EnuERecAvBias_untuned.Load(fin, kNOvAND, false, "EnuERecAvBias",
                                         kNuMu, "CH", selection(sel));

      auto splits_BANFF = Split(CH_numu_EnuERecAvBias.Target.get(), false);
      auto splits_NOvA = Split(CH_numu_EnuERecAvBias.From.get(), false);

      CH_numu_EnuERecAvBias_untuned.Target->SetName(
          "CH_numu_EnuERecAvBias_untuned_Target");
      CH_numu_EnuERecAvBias_untuned.From->SetName(
          "CH_numu_EnuERecAvBias_untuned_From");

      auto splits_BANFF_untuned =
          Split(CH_numu_EnuERecAvBias_untuned.Target.get(), false);
      auto splits_NOvA_untuned =
          Split(CH_numu_EnuERecAvBias_untuned.From.get(), false);

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

        splits_BANFF[i]->Scale(1.0 / GetMaximumBin(splits_BANFF[i]));
        splits_NOvA[i]->Scale(1.0 / GetMaximumBin(splits_NOvA[i]));
        splits_BANFF_untuned[i]->Scale(1.0 /
                                       GetMaximumBin(splits_BANFF_untuned[i]));
        splits_NOvA_untuned[i]->Scale(1.0 /
                                      GetMaximumBin(splits_NOvA_untuned[i]));

        splits_BANFF[i]->GetYaxis()->SetRangeUser(0, 1.01);

        splits_BANFF[i]->GetXaxis()->SetNdivisions(505);
        splits_BANFF[i]->GetYaxis()->SetNdivisions(505);

        splits_BANFF[i]->GetYaxis()->SetTitle("Relative rate");
        splits_BANFF[i]->GetYaxis()->SetTitleSize(0.06);
        splits_BANFF[i]->GetYaxis()->SetTitleOffset(0.8);
        splits_BANFF[i]->GetYaxis()->SetLabelSize(0.06);

        splits_BANFF[i]->GetXaxis()->SetTitleSize(0.06);
        splits_BANFF[i]->GetXaxis()->SetTitleOffset(1);
        splits_BANFF[i]->GetXaxis()->SetLabelSize(0.06);

        splits_BANFF[i]->SetLineColor(colorwheel[0]);
        splits_BANFF[i]->SetLineWidth(2);
        splits_BANFF[i]->SetLineStyle(1);
        splits_BANFF[i]->SetMarkerSize(0);
        splits_BANFF[i]->SetMarkerStyle(0);
        splits_BANFF[i]->SetMarkerColorAlpha(colorwheel[0], 0);

        splits_BANFF_untuned[i]->SetLineColor(colorwheel[1]);
        splits_BANFF_untuned[i]->SetLineWidth(2);
        splits_BANFF_untuned[i]->SetLineStyle(2);
        splits_BANFF_untuned[i]->SetMarkerSize(0);
        splits_BANFF_untuned[i]->SetMarkerStyle(0);
        splits_BANFF_untuned[i]->SetMarkerColorAlpha(colorwheel[1], 0);

        splits_NOvA[i]->SetLineColor(colorwheel[2]);
        splits_NOvA[i]->SetLineWidth(2);
        splits_NOvA[i]->SetLineStyle(1);
        splits_NOvA[i]->SetMarkerSize(0);
        splits_NOvA[i]->SetMarkerStyle(0);
        splits_NOvA[i]->SetMarkerColorAlpha(colorwheel[2], 0);

        splits_NOvA_untuned[i]->SetLineColor(colorwheel[3]);
        splits_NOvA_untuned[i]->SetLineWidth(2);
        splits_NOvA_untuned[i]->SetLineStyle(2);
        splits_NOvA_untuned[i]->SetMarkerSize(0);
        splits_NOvA_untuned[i]->SetMarkerStyle(0);
        splits_NOvA_untuned[i]->SetMarkerColorAlpha(colorwheel[1], 0);

        splits_BANFF[i]->DrawClone("HIST");
        splits_NOvA[i]->DrawClone("HISTSAME");
        splits_BANFF_untuned[i]->DrawClone("HISTSAME");
        splits_NOvA_untuned[i]->DrawClone("HISTSAME");

        c1.cd();
      }

      TLegend *leg = nullptr;

      if (ybuffer > 0) {
        leg = new TLegend(0.1, 0.85, 0.9, 1);
        leg->SetTextSize(0.04);
      } else { // Put the legend in the place of the last pane(s)
        leg = new TLegend((nslices % nx) * padxwidth, 1 - (ny * padywidth),
                          nx * padxwidth, 1 - ((ny - 1) * padywidth));
        leg->SetTextSize(0.04 * float((nx * ny) - nslices));
      }

      leg->SetBorderSize(0);
      leg->SetFillStyle(0);
      leg->SetNColumns(2);
      leg->AddEntry(splits_BANFF[first_slice], "BANFF", "lp");
      leg->AddEntry(splits_BANFF_untuned[first_slice], "NEUT", "lp");
      leg->AddEntry(splits_NOvA[first_slice], "NOvA", "lp");
      leg->AddEntry(splits_NOvA_untuned[first_slice], "GENIE", "lp");

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