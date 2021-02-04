
#include "T2KNOvAFakeDataHelper.hxx"

struct hblob {

  TH1 *From;
  TH1 *Target;
  std::vector<TH1 *> ReWeights;

  bool ist2kbase;

  void Load(TFile *fin, bool t2kbase, std::string varname,
            t2knova::nuspecies nuspec, t2knova::selection sel) {
    this->ist2kbase = t2kbase;

    From = GetTH1(fin, std::string(t2kbase ? "NEUT/ND280/" : "GENIE/NOvAND/") +
                           t2knova::all_nuspecies[nuspec] + "/norw/" + varname +
                           "_" + t2knova::all_sel[sel]);

    Target =
        GetTH1(fin, std::string(t2kbase ? "GENIE/ND280/" : "NEUT/NOvAND/") +
                        t2knova::all_nuspecies[nuspec] + "/norw/" + varname +
                        "_" + t2knova::all_sel[sel]);

    if (t2kbase) {
      ReWeights.push_back(GetTH1(
          fin, std::string("NEUT/ND280/") + t2knova::all_nuspecies[nuspec] +
                   "/rw_GENIE/" + varname + "_" + t2knova::all_sel[sel]));
    } else {
      ReWeights.push_back(GetTH1(
          fin, std::string("GENIE/NOvAND/") + t2knova::all_nuspecies[nuspec] +
                   "/rw_plep_NEUT/" + varname + "_" + t2knova::all_sel[sel]));
      ReWeights.push_back(GetTH1(
          fin, std::string("GENIE/NOvAND/") + t2knova::all_nuspecies[nuspec] +
                   "/rw_Q2_NEUT/" + varname + "_" + t2knova::all_sel[sel]));
      ReWeights.push_back(GetTH1(
          fin, std::string("GENIE/NOvAND/") + t2knova::all_nuspecies[nuspec] +
                   "/rw_ptlep_NEUT/" + varname + "_" + t2knova::all_sel[sel]));
    }
  }

  void Print(const char *fname, const char *title = "") {

    double max_gen = GetMaximumTH1s({From, Target});
    double max_rw = GetMaximumTH1s(ReWeights);
    double max = std::max(max_gen, max_rw);

    TCanvas *c1 = MakeCanvasTopLegend();

    TPad *p1 = MakeRatioTopPadTopLegend();
    p1->AppendPad();
    TPad *p2 = MakeRatioBottomPadTopLegend();
    p2->AppendPad();

    p1->cd();

    From->SetLineColor(kBlack);
    From->SetLineWidth(2);

    Target->SetLineWidth(2);
    Target->SetLineColor(SORNBrightWheel[1]);

    int cols[] = {SORNBrightWheel[2], SORNBrightWheel[4], SORNBrightWheel[5]};
    int lc = 0;
    for (auto &h : ReWeights) {
      h->SetLineWidth(2);
      h->SetLineColor(cols[lc++]);
    }

    From->GetYaxis()->SetRangeUser(0, max * 1.1);
    From->GetYaxis()->SetTitle("Cross Section 10^{-39}");
    StyleAxis(From->GetYaxis());
    HideAxis(From->GetXaxis());

    From->DrawClone("EHIST");
    Target->DrawClone("EHISTSAME");

    for (auto &h : ReWeights) {
      h->DrawClone("EHISTSAME");
    }

    p2->cd();

    From->Divide(Target);
    for (auto &h : ReWeights) {
      h->Divide(Target);
    }

    Target->Divide(Target);
    StyleAxis(Target->GetXaxis(), 2);
    StyleAxis(Target->GetYaxis(), 2, 0.5, 1);
    Target->GetYaxis()->SetRangeUser(0.8, 1.2);
    Target->GetYaxis()->SetTitle("Ratio To Generated");
    Target->SetLineStyle(2);

    Target->Draw("HIST");
    From->DrawClone("HISTSAME");

    for (auto &h : ReWeights) {
      h->DrawClone("HISTSAME");
    }

    c1->cd();

    Target->SetLineStyle(1);

    TLegend *leg = MakeTopLegend();
    leg->SetTextSize(0.03);

    leg->AddEntry(From, ist2kbase ? "BANFF Post ND280" : "NOvA2020", "l");
    leg->AddEntry(Target, ist2kbase ? "NOvA2020" : "BANFF Post ND280", "l");
    leg->AddEntry(ReWeights[0],
                  ist2kbase ? "Reweight to NOvA" : "ReWeight to BANFF (PLep)",
                  "l");
    if (!ist2kbase) {
      leg->AddEntry(ReWeights[1], "ReWeight to BANFF (Q^{2})", "l");
      leg->AddEntry(ReWeights[2], "ReWeight to BANFF (PtLep)", "l");
    }

    leg->Draw();

    TLatex ttl;
    ttl.SetTextSize(0.05);
    ttl.SetTextAlign(22);
    ttl.DrawLatexNDC(0.75, 0.825, title);

    c1->Print(fname);
  }

  static void LoadAndPrint(TFile *fin, bool t2kbase, std::string varname,
                           t2knova::nuspecies nuspec, t2knova::selection sel,
                           const char *fname, const char *title) {
    hblob h;
    h.Load(fin, t2kbase, varname, nuspec, sel);
    h.Print(fname, title);
  }
};

void ValidPlots(std::string const &finname) {

  TFile fin(finname.c_str());
  if (fin.IsZombie()) {
    std::cout << "Failed to read " << finname << std::endl;
    return;
  }

  TCanvas c1("cdummy", "", 800, 800);
  c1.Print("validplots.pdf[");

  hblob::LoadAndPrint(&fin, true, "Enu", t2knova::kNuMu, t2knova::kCCINC,
                      "validplots.pdf", "CCInc");
  hblob::LoadAndPrint(&fin, true, "Enu", t2knova::kNuMu, t2knova::kCC0pi,
                      "validplots.pdf", "CC0#pi");
  hblob::LoadAndPrint(&fin, true, "Enu", t2knova::kNuMu, t2knova::kCC1cpi,
                      "validplots.pdf", "CC1#pi^{#pm}");
  hblob::LoadAndPrint(&fin, true, "Enu", t2knova::kNuMu, t2knova::kCC1pi0,
                      "validplots.pdf", "CC1#pi^{0}");
    hblob::LoadAndPrint(&fin, true, "Enu", t2knova::kNuMu, t2knova::kCCOther,
                      "validplots.pdf", "CCOther");

  hblob::LoadAndPrint(&fin, false, "Enu", t2knova::kNuMu, t2knova::kCCINC,
                      "validplots.pdf", "CCInc");
  hblob::LoadAndPrint(&fin, false, "Enu", t2knova::kNuMu, t2knova::kCC0pi,
                      "validplots.pdf", "CC0#pi");
  hblob::LoadAndPrint(&fin, false, "Enu", t2knova::kNuMu, t2knova::kCC1cpi,
                      "validplots.pdf", "CC1#pi^{#pm}");
  hblob::LoadAndPrint(&fin, false, "Enu", t2knova::kNuMu, t2knova::kCC1pi0,
                      "validplots.pdf", "CC1#pi^{0}");
  hblob::LoadAndPrint(&fin, false, "Enu", t2knova::kNuMu, t2knova::kCCOther,
                      "validplots.pdf", "CCOther");

  hblob::LoadAndPrint(&fin, true, "PLep", t2knova::kNuMu, t2knova::kCCINC,
                      "validplots.pdf", "CCInc");
  hblob::LoadAndPrint(&fin, true, "PLep", t2knova::kNuMu, t2knova::kCC0pi,
                      "validplots.pdf", "CC0#pi");
  hblob::LoadAndPrint(&fin, true, "PLep", t2knova::kNuMu, t2knova::kCC1cpi,
                      "validplots.pdf", "CC1#pi^{#pm}");
  hblob::LoadAndPrint(&fin, true, "PLep", t2knova::kNuMu, t2knova::kCC1pi0,
                      "validplots.pdf", "CC1#pi^{0}");
    hblob::LoadAndPrint(&fin, true, "PLep", t2knova::kNuMu, t2knova::kCCOther,
                      "validplots.pdf", "CCOther");

  hblob::LoadAndPrint(&fin, false, "PLep", t2knova::kNuMu, t2knova::kCCINC,
                      "validplots.pdf", "CCInc");
  hblob::LoadAndPrint(&fin, false, "PLep", t2knova::kNuMu, t2knova::kCC0pi,
                      "validplots.pdf", "CC0#pi");
  hblob::LoadAndPrint(&fin, false, "PLep", t2knova::kNuMu, t2knova::kCC1cpi,
                      "validplots.pdf", "CC1#pi^{#pm}");
  hblob::LoadAndPrint(&fin, false, "PLep", t2knova::kNuMu, t2knova::kCC1pi0,
                      "validplots.pdf", "CC1#pi^{0}");
  hblob::LoadAndPrint(&fin, false, "PLep", t2knova::kNuMu, t2knova::kCCOther,
                      "validplots.pdf", "CCOther");

  hblob::LoadAndPrint(&fin, true, "Q2", t2knova::kNuMu, t2knova::kCCINC,
                      "validplots.pdf", "CCInc");
  hblob::LoadAndPrint(&fin, true, "Q2", t2knova::kNuMu, t2knova::kCC0pi,
                      "validplots.pdf", "CC0#pi");
  hblob::LoadAndPrint(&fin, true, "Q2", t2knova::kNuMu, t2knova::kCC1cpi,
                      "validplots.pdf", "CC1#pi^{#pm}");
  hblob::LoadAndPrint(&fin, true, "Q2", t2knova::kNuMu, t2knova::kCC1pi0,
                      "validplots.pdf", "CC1#pi^{0}");
    hblob::LoadAndPrint(&fin, true, "Q2", t2knova::kNuMu, t2knova::kCCOther,
                      "validplots.pdf", "CCOther");

  hblob::LoadAndPrint(&fin, false, "Q2", t2knova::kNuMu, t2knova::kCCINC,
                      "validplots.pdf", "CCInc");
  hblob::LoadAndPrint(&fin, false, "Q2", t2knova::kNuMu, t2knova::kCC0pi,
                      "validplots.pdf", "CC0#pi");
  hblob::LoadAndPrint(&fin, false, "Q2", t2knova::kNuMu, t2knova::kCC1cpi,
                      "validplots.pdf", "CC1#pi^{#pm}");
  hblob::LoadAndPrint(&fin, false, "Q2", t2knova::kNuMu, t2knova::kCC1pi0,
                      "validplots.pdf", "CC1#pi^{0}");
  hblob::LoadAndPrint(&fin, false, "Q2", t2knova::kNuMu, t2knova::kCCOther,
                      "validplots.pdf", "CCOther");

  hblob::LoadAndPrint(&fin, true, "PtLep", t2knova::kNuMu, t2knova::kCCINC,
                      "validplots.pdf", "CCInc");
  hblob::LoadAndPrint(&fin, true, "PtLep", t2knova::kNuMu, t2knova::kCC0pi,
                      "validplots.pdf", "CC0#pi");
  hblob::LoadAndPrint(&fin, true, "PtLep", t2knova::kNuMu, t2knova::kCC1cpi,
                      "validplots.pdf", "CC1#pi^{#pm}");
  hblob::LoadAndPrint(&fin, true, "PtLep", t2knova::kNuMu, t2knova::kCC1pi0,
                      "validplots.pdf", "CC1#pi^{0}");
    hblob::LoadAndPrint(&fin, true, "PtLep", t2knova::kNuMu, t2knova::kCCOther,
                      "validplots.pdf", "CCOther");

  hblob::LoadAndPrint(&fin, false, "PtLep", t2knova::kNuMu, t2knova::kCCINC,
                      "validplots.pdf", "CCInc");
  hblob::LoadAndPrint(&fin, false, "PtLep", t2knova::kNuMu, t2knova::kCC0pi,
                      "validplots.pdf", "CC0#pi");
  hblob::LoadAndPrint(&fin, false, "PtLep", t2knova::kNuMu, t2knova::kCC1cpi,
                      "validplots.pdf", "CC1#pi^{#pm}");
  hblob::LoadAndPrint(&fin, false, "PtLep", t2knova::kNuMu, t2knova::kCC1pi0,
                      "validplots.pdf", "CC1#pi^{0}");
  hblob::LoadAndPrint(&fin, false, "PtLep", t2knova::kNuMu, t2knova::kCCOther,
                      "validplots.pdf", "CCOther");

  c1.Print("validplots.pdf]");
}