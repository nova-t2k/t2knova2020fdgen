#include "T2KNOvAFakeDataHelper.hxx"
#include "TLatex.h"
#include "colordef.h"
#include "plotutils.h"

struct hblob {
  std::unique_ptr<TH1> From;
  std::unique_ptr<TH1> Target;
  std::vector<std::unique_ptr<TH1>> ReWeights;

  bool ist2kbase;

  void Load(std::unique_ptr<TFile> &fin, bool t2kbase, std::string varname,
            t2knova::nuspecies nuspec, std::string tgtstr,
            t2knova::selection sel) {
    this->ist2kbase = t2kbase;

    From = GetTH1(
        fin, std::string(t2kbase ? "ND280/T2KNDTune" : "NOvAND/NOvATune") +
                 "/" + tgtstr + "/" + t2knova::all_nuspecies[nuspec] + "/" +
                 varname + "_" + t2knova::SelectionList[sel]);

    Target = GetTH1(
        fin, std::string(t2kbase ? "ND280/NOvATune" : "NOvAND/T2KNDTune") +
                 "/" + tgtstr + "/" + t2knova::all_nuspecies[nuspec] + "/" +
                 varname + "_" + t2knova::SelectionList[sel]);

    if (t2kbase) {
      ReWeights.push_back(
          GetTH1(fin, std::string("ND280/T2KNDTune_To_NOvATune") + "/" +
                          tgtstr + "/" + t2knova::all_nuspecies[nuspec] + "/" +
                          varname + "_" + t2knova::SelectionList[sel]));
      ReWeights.push_back(
          GetTH1(fin, std::string("ND280/T2KNDTune_To_NOvATune_Enu") + "/" +
                          tgtstr + "/" + t2knova::all_nuspecies[nuspec] + "/" +
                          varname + "_" + t2knova::SelectionList[sel]));
      ReWeights.push_back(
          GetTH1(fin, std::string("ND280/T2KNDTune_To_NOvATune_Q2") + "/" +
                          tgtstr + "/" + t2knova::all_nuspecies[nuspec] + "/" +
                          varname + "_" + t2knova::SelectionList[sel]));
    } else {
      throw;
    }
  }

  void Print(const char *fname, const char *title = "") {
    if (!From || !Target) {
      return;
    }

    double max_gen = GetMaximumTH1s(std::vector<std::reference_wrapper<std::unique_ptr<TH1>>>{From, Target});
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

    From->Divide(Target.get());
    for (auto &h : ReWeights) {
      h->Divide(Target.get());
    }

    Target->Divide(Target.get());
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

    leg->AddEntry(From.get(), ist2kbase ? "BANFF Post ND280" : "NOvA2020", "l");
    leg->AddEntry(Target.get(), ist2kbase ? "NOvA2020" : "BANFF Post ND280",
                  "l");
    leg->AddEntry(ReWeights[0].get(),
                  ist2kbase ? "Reweight to NOvA (EnuPThetaLep)"
                            : "ReWeight to BANFF (PLep)",
                  "l");
    if (!ist2kbase) {
      leg->AddEntry(ReWeights[1].get(), "ReWeight to BANFF (Q^{2})", "l");
      leg->AddEntry(ReWeights[2].get(), "ReWeight to BANFF (PtLep)", "l");
    } else {
      leg->AddEntry(ReWeights[1].get(), "ReWeight to NOvA (Enu)", "l");
      leg->AddEntry(ReWeights[2].get(), "ReWeight to NOvA (Q^{2})", "l");
    }

    leg->Draw();

    TLatex ttl;
    ttl.SetTextSize(0.05);
    ttl.SetTextAlign(22);
    ttl.DrawLatexNDC(0.75, 0.825, title);

    c1->Print(fname);
  }

  static void LoadAndPrint(std::unique_ptr<TFile> &fin, bool t2kbase,
                           std::string varname, t2knova::nuspecies nuspec,
                           std::string tgtstr, t2knova::selection sel,
                           const char *fname, const char *title) {
    hblob h;
    h.Load(fin, t2kbase, varname, nuspec, tgtstr, sel);
    h.Print(fname, (std::string(title) + " " + tgtstr + " " +
                    t2knova::all_nuspecies_latex[nuspec])
                       .c_str());
  }
};

void ValidPlots(std::string const &finname) {
  std::unique_ptr<TFile> fin(new TFile(finname.c_str()));
  if (fin->IsZombie()) {
    std::cout << "Failed to read " << finname << std::endl;
    return;
  }

  TCanvas c1("cdummy", "", 800, 800);
  c1.Print("validplots.pdf[");

  for (auto tgtstr : {"C", "H", "O", "CH", "H2O"}) {
    for (auto nuspec :
         {t2knova::kNuMu, t2knova::kNuMub, t2knova::kNuE, t2knova::kNuEb}) {
      for (auto proj : {
               "Enu",
               "PLep",
               "Q2",
               "PtLep",
           })
        for (int sel = 0; sel < t2knova::SelectionList.size(); ++sel) {
          hblob::LoadAndPrint(fin, true, proj, nuspec, tgtstr,
                              t2knova::selection(sel), "validplots.pdf",
                              t2knova::SelectionList[sel].c_str());
        }
    }
  }

  c1.Print("validplots.pdf]");
}

int main(int argc, char const *argv[]) {
  DeclareColors();
  gStyle->SetOptStat(false);
  ValidPlots(argv[1]);
}