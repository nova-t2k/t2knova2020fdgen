#include "T2KNOvA/FakeDataHelper.hxx"
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
            t2knova::selection sel, int mode = 0) {
    this->ist2kbase = t2kbase;

    std::string mode_str = "";
    if (mode != 0) {
      mode_str =
          (mode < 0 ? "_Mode_m" : "_Mode_") + std::to_string(std::abs(mode));
    }

    From =
        GetTH1(fin,
               std::string(t2kbase ? "ND280/T2KNDTune" : "NOvAND/NOvATune") +
                   "/" + tgtstr + "/" + t2knova::all_nuspecies[nuspec] + "/" +
                   varname + "_" + t2knova::SelectionList[sel] + mode_str,
               false);

    Target =
        GetTH1(fin,
               std::string(t2kbase ? "ND280/NOvATune" : "NOvAND/T2KNDTune") +
                   "/" + tgtstr + "/" + t2knova::all_nuspecies[nuspec] + "/" +
                   varname + "_" + t2knova::SelectionList[sel] + mode_str,
               false);

    if (!From || !Target) {
      return;
    }

    From->SetName("FROM");
    Target->SetName("TARGET");

    // std::vector<std::string> ReWeightHists = ist2kbase ?
    // std::vector<std::string>{
    //   std::string("ND280/T2KND_to_NOvA") + "/" + tgtstr + "/" +
    //                  t2knova::all_nuspecies[nuspec] + "/" + varname + "_" +
    //                  t2knova::SelectionList[sel],
    //                  std::string("ND280/T2KND_to_NOvA_EnuKludge") + "/" +
    //                  tgtstr +
    //                  "/" + t2knova::all_nuspecies[nuspec] + "/" + varname +
    //                  "_" + t2knova::SelectionList[sel],
    //                  std::string("ND280/T2KND_to_NOvA_Enu") + "/" + tgtstr +
    //                  "/" + t2knova::all_nuspecies[nuspec] + "/" + varname +
    //                  "_" + t2knova::SelectionList[sel],
    //                  std::string("ND280/T2KND_to_NOvA_Q2") + "/" + tgtstr +
    //                  "/" + t2knova::all_nuspecies[nuspec] + "/" + varname +
    //                  "_" + t2knova::SelectionList[sel],
    // } : std::vector<std::string>{
    //     std::string("NOvAND/NOvA_to_T2KND_ptlep") + "/" + tgtstr + "/" +
    //                  t2knova::all_nuspecies[nuspec] + "/" + varname + "_" +
    //                  t2knova::SelectionList[sel],
    //                  std::string("NOvAND/NOvA_to_T2KND_plep") + "/" + tgtstr
    //                  + "/" + t2knova::all_nuspecies[nuspec] + "/" + varname +
    //                  "_" + t2knova::SelectionList[sel],
    //                  std::string("NOvAND/NOvA_to_T2KND_Q2") + "/" + tgtstr +
    //                  "/" + t2knova::all_nuspecies[nuspec] + "/" + varname +
    //                  "_" + t2knova::SelectionList[sel],
    // };

    std::vector<std::string> ReWeightHists = ist2kbase ? std::vector<std::string>{
      std::string("ND280/T2KND_to_NOvA") + "/" + tgtstr + "/" +
                     t2knova::all_nuspecies[nuspec] + "/" + varname + "_" +
                     t2knova::SelectionList[sel] + mode_str,} : std::vector<std::string>{
        std::string("NOvAND/NOvA_to_T2KND_ptlep") + "/" + tgtstr + "/" +
                     t2knova::all_nuspecies[nuspec] + "/" + varname + "_" +
                     t2knova::SelectionList[sel] + mode_str,};

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

  void Print(const char *fname, const char *title = "") {
    if (!From || !Target) {
      return;
    }

    double max_gen = GetMaximumTH1s(
        std::vector<std::reference_wrapper<std::unique_ptr<TH1>>>{From,
                                                                  Target});
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

    From->GetYaxis()->SetRangeUser(0, max * 1.1);
    From->GetYaxis()->SetTitle("Cross Section 10^{-39}");
    From->SetTitle("");
    StyleAxis(From->GetYaxis());
    HideAxis(From->GetXaxis());

    ReWeights[0]->Draw();

    From->DrawClone("EHIST");
    Target->DrawClone("EHISTSAME");

    for (auto &h : ReWeights) {
      if (!h) {
        throw;
        continue;
      }
      h->DrawClone("EHISTSAME");
    }

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

    c1->cd();

    Target->SetLineStyle(1);

    TLegend *leg = MakeTopLegend();
    leg->SetTextSize(0.02);

    leg->AddEntry(From.get(), ist2kbase ? "BANFF Post ND280" : "NOvA2020", "l");
    leg->AddEntry(Target.get(), ist2kbase ? "NOvA2020" : "BANFF Post ND280",
                  "l");

    if (ReWeights[0]) {
      leg->AddEntry(ReWeights[0].get(),
                    ist2kbase ? "Reweight to NOvA (Enu PLep ThetaLep)"
                              : "ReWeight to BANFF (Enu PtLep EVisHad)",
                    "l");
    }

    leg->Draw();

    TLatex ttl;
    ttl.SetTextSize(0.05);
    ttl.SetTextAlign(22);
    ttl.DrawLatexNDC(0.75, 0.75, title);

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

    bool bymode = true;
    if (bymode) {
      for (int i = -60; i < 60; ++i) {
        if(i == 0){
          continue;
        }
        hblob h;
        h.Load(fin, t2kbase, varname, nuspec, tgtstr, sel, i);
        h.Print(fname, (std::string(title) + " " + tgtstr + " " +
                        t2knova::all_nuspecies_latex[nuspec] +
                        " m:" + std::to_string(i))
                           .c_str());
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
  c1.Print("validplots_t2k.pdf[");
  c1.Print("validplots_nova.pdf[");

  // for (auto tgtstr : {"C", "H", "O", "CH", "H2O"}) {
  for (auto tgtstr : {std::string("CH"), std::string("H2O")}) {
    for (auto nuspec :
         {t2knova::kNuMu, t2knova::kNuMub, t2knova::kNuE, t2knova::kNuEb}) {
      for (int sel : {
               t2knova::kCCInc,
               t2knova::kCC0pi,
               t2knova::kCC1cpi,
               t2knova::kCC1pi0,
               t2knova::kCCmultipi,
               t2knova::kCC1Gamma,
               t2knova::kCCOther,
               t2knova::kNCInc,
               t2knova::kNC0pi,
               t2knova::kNC1cpi,
               t2knova::kNC1pi0,
               t2knova::kNCmultipi,
               t2knova::kNC1Gamma,
               t2knova::kNCOther,
           }) {
        for (auto proj : {
                 "Enu",
                 "ERecQE",
                 "PLep",
                 "ThetaLep",
                 "Q2",
                 "q0",
                 "q3",
                 "hmfscpip",
                 "hmfspi0p",
                 "ncpi",
                 "npi0",
             }) {

          if(std::string(proj) == "ERecQE" && sel >= t2knova::kNCInc ){
            continue;
          }
          hblob::LoadAndPrint(fin, true, proj, nuspec, tgtstr,
                              t2knova::selection(sel), "validplots_t2k.pdf",
                              t2knova::SelectionList[sel].c_str());
        }
        // for (auto proj : {
        //          "Enu",
        //          "PtLep",
        //          "PLep",
        //          "EAvHad",
        //          "Q2",
        //          "q0",
        //          "q3",
        //          "hmfscpip",
        //          "hmfspi0p",
        //          "ncpi",
        //          "npi0",
        //      }) {
        //   hblob::LoadAndPrint(fin, false, proj, nuspec, tgtstr,
        //                       t2knova::selection(sel), "validplots_nova.pdf",
        //                       t2knova::SelectionList[sel].c_str());
        // }
      }
    }
  }

  c1.Print("validplots_t2k.pdf]");
  c1.Print("validplots_nova.pdf]");
}

int main(int argc, char const *argv[]) {
  DeclareColors();
  gStyle->SetOptStat(false);
  ValidPlots(argv[1]);
}