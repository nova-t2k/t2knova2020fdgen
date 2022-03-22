#include "T2KNOvA/FakeDataHelper.hxx"
#include "TLatex.h"
#include "colordef.h"
#include "plotutils.h"

#include "TExec.h"

using namespace t2knova;

bool const kT2K = true;
bool const kNOvA = false;
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

struct hblob {
  std::unique_ptr<TH1> From;
  std::unique_ptr<TH1> Target;
  std::vector<std::unique_ptr<TH1>> ReWeights;

  bool is2D;
  bool ist2kbase;

  void Load(std::unique_ptr<TFile> &fin, bool t2kbase, bool isosc,
            std::string varname, nuspecies nuspec, std::string tgtstr,
            selection sel, int mode = 0) {
    this->ist2kbase = t2kbase;

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

    From =
        GetTH1(fin,
               std::string(t2kbase ? "ND280/T2KNDTune" : "NOvAND/NOvATune") +
                   (isosc ? "_osc" : "") + "/" + tgtstr + "/" +
                   all_nuspecies[nuspec] + "/" + varname + sel_str + mode_str,
               false);

    Target =
        GetTH1(fin,
               std::string(t2kbase ? "ND280/NOvATune" : "NOvAND/T2KNDTune") +
                   (isosc ? "_osc" : "") + "/" + tgtstr + "/" +
                   all_nuspecies[nuspec] + "/" + varname + sel_str + mode_str,
               false);

    if (!From || !Target) {
      return;
    }

    From->SetName("FROM");
    Target->SetName("TARGET");

    std::vector<std::string> ReWeightHists;

    if (ist2kbase) {
      ReWeightHists = std::vector<std::string>{
          std::string("ND280/T2KND_to_NOvA") + (isosc ? "_osc" : "") + "/" +
              tgtstr + "/" + all_nuspecies[nuspec] + "/" + varname + sel_str +
              mode_str,
      };
    } else {
      ReWeightHists = std::vector<std::string>{
          std::string("NOvAND/NOvA_to_T2KND_ptlep") + (isosc ? "_osc" : "") +
              "/" + tgtstr + "/" + all_nuspecies[nuspec] + "/" + varname +
              sel_str + mode_str,
      };
    }

    int i = 0;
    for (auto &n : ReWeightHists) {
      ReWeights.push_back(GetTH1(fin, n));
      if (!ReWeights.back()) {
        std::cout << "[ERROR]: Failed to read " << n << std::endl;
      }
      ReWeights.back()->SetName(
          (std::string("ReWeight_") + std::to_string(i++)).c_str());
      std::cout << "read " << n << " as " << ReWeights.back()->GetName()
                << std::endl;
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
      From->SetLineWidth(2);
    }

    if (Target) {
      Target->SetLineWidth(2);
      Target->SetLineColor(SORNBrightWheel[1]);
    }

    int cols[] = {SORNBrightWheel[2], SORNBrightWheel[4], SORNBrightWheel[5],
                  SORNBrightWheel[6]};
    int lc = 0;
    for (auto &h : ReWeights) {
      if (!h) {
        lc++;
        continue;
      }
      h->SetLineWidth(2);
      h->SetLineStyle(2);
      h->SetLineColor(cols[lc++]);
    }

    First->GetYaxis()->SetRangeUser(0, max * 1.1);
    First->GetYaxis()->SetTitle("Cross Section 10^{-39}");
    First->SetTitle("");
    StyleAxis(First->GetYaxis());

    if (HaveBoth) {
      HideAxis(First->GetXaxis());
    }

    if (HaveReWeight) {
      ReWeights[0]->Draw();
    }

    First->DrawClone("EHIST");
    if (HaveBoth) {
      Target->DrawClone("EHISTSAME");
    }

    for (auto &h : ReWeights) {
      if (!h) {
        throw;
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

      Target->Divide(Target.get());
      StyleAxis(Target->GetXaxis(), 2);
      StyleAxis(Target->GetYaxis(), 2, 0.5, 1);
      Target->GetYaxis()->SetRangeUser(0.8, 1.2);
      Target->GetYaxis()->SetTitle("Ratio To Generated");
      Target->SetTitle("");
      Target->SetLineStyle(2);

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

    TLegend *leg = MakeTopLegend();
    leg->SetTextSize(0.02);

    leg->AddEntry(From.get(), ist2kbase ? "BANFF Post ND280" : "NOvA2020", "l");
    leg->AddEntry(Target.get(), ist2kbase ? "NOvA2020" : "BANFF Post ND280",
                  "l");

    if (HaveReWeight) {
      leg->AddEntry(ReWeights[0].get(),
                    ist2kbase ? "Reweight to NOvA (Enu PLep ThetaLep)"
                              : "ReWeight to BANFF (Enu PtLep EVisHad)",
                    "l");
    }

    leg->Draw();

    TLatex ttl;
    ttl.SetTextSize(0.04);
    ttl.SetTextAlign(32);
    ttl.DrawLatexNDC(0.99, 0.875, title.c_str());
    if (mode_line.size()) {
      ttl.DrawLatexNDC(0.99, 0.825, mode_line.c_str());
    }

    OpenPDF(fname);
    c1->Print((fname + ".pdf").c_str());
  }

  void Print2D(std::string fname, std::string title = "",
               std::string mode_line = "") {

    static const long unsigned int BWRPalette =
        (long unsigned int)GetBWRPalette();

    static std::string BWR_cmd = "gStyle->SetPalette(100, (int*)" +
                                 std::to_string(BWRPalette) +
                                 ");gStyle->SetNumberContours(100);";

    static TExec *ex1 = new TExec("ex1", "gStyle->SetPalette(kBird);");
    static TExec *ex2 = new TExec("ex2", BWR_cmd.c_str());

    bool HaveBoth = bool(From) && bool(Target);
    bool HaveReWeight = ReWeights.size() && bool(ReWeights[0]);

    TCanvas *c1 = MakeCanvas();

    TLatex ttl;
    ttl.SetTextSize(0.025);
    ttl.SetTextAlign(12);

    if (From) {
      TPad *p = new TPad("pFrom", "", 0, 0.5, 0.5, 1);
      p->AppendPad();
      p->SetTopMargin(0.05);
      p->SetBottomMargin(0.2);
      p->SetLeftMargin(0.11);
      p->SetRightMargin(0.14);
      p->SetTopMargin(0.03);
      p->cd();

      From->SetTitle("");
      auto CL = From->DrawClone("COLZ");
      ex1->Draw();
      CL->Draw("COLZ SAME");

      c1->cd();
      ttl.DrawLatexNDC(0.1, 0.95, ist2kbase ? "BANFF Post ND280" : "NOvA2020");
    }

    if (Target) {
      TPad *p = new TPad("pTarget", "", 0.5, 0.5, 1, 1);
      p->AppendPad();
      p->SetTopMargin(0.05);
      p->SetBottomMargin(0.2);
      p->SetLeftMargin(0.11);
      p->SetRightMargin(0.14);
      p->SetTopMargin(0.03);
      p->cd();

      Target->SetTitle("");
      auto CL = Target->DrawClone("COLZ");
      ex1->Draw();
      CL->Draw("COLZ SAME");

      c1->cd();
      ttl.DrawLatexNDC(0.6, 0.95, ist2kbase ? "NOvA2020" : "BANFF Post ND280");
    }

    if (HaveBoth) {
      From->Divide(Target.get());

      TPad *p = new TPad("pBoth", "", 0, 0, 0.5, 0.5);
      p->AppendPad();
      p->SetTopMargin(0.05);
      p->SetBottomMargin(0.2);
      p->SetLeftMargin(0.11);
      p->SetRightMargin(0.14);
      p->SetTopMargin(0.03);
      p->cd();

      From->GetZaxis()->SetRangeUser(0, 2);
      From->GetZaxis()->SetTitleOffset(1.25);
      From->GetZaxis()->SetTitle("Event Rate Ratio");
      auto CL = From->DrawClone("COLZ");
      ex2->Draw();
      CL->Draw("COLZ SAME");

      c1->cd();
      ttl.DrawLatexNDC(0.1, 0.45,
                       ist2kbase ? "BANFF/NOvA2020" : "NOvA2020/ND280");
    }

    if (HaveReWeight) {
      TPad *p = new TPad("pRW", "", 0.5, 0, 1, 0.5);
      p->AppendPad();
      p->SetTopMargin(0.05);
      p->SetBottomMargin(0.2);
      p->SetLeftMargin(0.11);
      p->SetRightMargin(0.14);
      p->SetTopMargin(0.03);
      p->cd();

      ReWeights[0]->Divide(Target.get());
      ReWeights[0]->SetTitle("");
      ReWeights[0]->GetZaxis()->SetRangeUser(0, 2);
      ReWeights[0]->GetZaxis()->SetTitleOffset(1.25);
      ReWeights[0]->GetZaxis()->SetTitle("Event Rate Ratio");
      auto CL = ReWeights[0]->DrawClone("COLZ");
      ex2->Draw();
      CL->Draw("COLZ SAME");

      c1->cd();
      ttl.DrawLatexNDC(0.6, 0.45, "Reweight/Generated");
    }
    OpenPDF(fname);
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

  static void LoadAndPrint(std::unique_ptr<TFile> &fin, bool t2kbase,
                           std::string varname, nuspecies nuspec,
                           std::string tgtstr, selection sel, std::string fname,
                           std::string title) {
    hblob h;
    h.Load(fin, t2kbase, false, varname, nuspec, tgtstr, sel);
    h.Print(fname, title + " " + tgtstr + " " + all_nuspecies_latex[nuspec],
            "All modes");

    if (doosc) {
      h.Load(fin, t2kbase, true, varname, nuspec, tgtstr, sel);
      h.Print(fname + "_osc",
              title + " " + tgtstr + " " + all_nuspecies_latex[nuspec],
              "All modes");
    }

    if (bymode) {
      for (int i = -60; i < 60; ++i) {
        if (i == 0) {
          continue;
        }
        h.Load(fin, t2kbase, false, varname, nuspec, tgtstr, sel, i);
        h.Print(fname + "_Modes",
                title + " " + tgtstr + " " + all_nuspecies_latex[nuspec],
                "True mode:" + std::to_string(i));

        if (doosc) {
          h.Load(fin, t2kbase, true, varname, nuspec, tgtstr, sel, i);
          h.Print(fname + "_Modes_osc",
                  title + " " + tgtstr + " " + all_nuspecies_latex[nuspec],
                  "True mode:" + std::to_string(i));
        }
      }
    }
  }
};

void ValidPlots(std::string const &finname) {
  std::unique_ptr<TFile> fin(new TFile(finname.c_str()));
  if (fin->IsZombie()) {
    std::cout << "Failed to read " << finname << std::endl;
    return;
  }

  TCanvas c1("cdummy", "", 800, 800);

  for (auto tgtstr : {std::string("CH"), std::string("H2O")}) {
    for (auto nuspec : {kNuMu, kNuMub, kNuE, kNuEb}) {

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

        for (auto t2k_proj : {"Enu", "ERecQE", "PLep", "ThetaLep", "Q2",
                              "q0q3_low", "q0q3_high", "hmfscpip", "hmfspi0p",
                              "ncpi", "npi0", "EGamma", "EGamma_DeExcite"}) {

          hblob::LoadAndPrint(
              fin, kT2K, t2k_proj, nuspec, tgtstr, selection(sel),
              "Valid_ND280_" + SelectionList[sel] + "_" + all_nuspecies[nuspec],
              SelectionList[sel]);
        }
      }
    }
  }

  ClosePDFs();
}

int main(int argc, char const *argv[]) {
  DeclareColors();
  gStyle->SetOptStat(false);
  ValidPlots(argv[1]);
}