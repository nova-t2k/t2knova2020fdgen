#include "TH2.h"
#include "TH3.h"
#include "TFile.h"

#include <string>
#include <iostream>
#include <cmath>

TH1 *GetTH1(TFile *f, std::string const &name) {
  TDirectory *odir = gDirectory;

  TH1 *h;
  f->GetObject(name.c_str(), h);
  if (!h) {
    std::cout << "[ERROR]: Failed to find histogram named: " << name
              << std::endl;
  } else {
    h->SetDirectory(nullptr);
  }

  if (odir) {
    gDirectory->cd();
  }

  return h;
}

void ScrubLowStatsBins(TH3D *num, TH3D *denom, TH3D *ratio,
                       double frac_error_threshold) {
  for (int i = 0; i < num->GetXaxis()->GetNbins(); ++i) {
    for (int j = 0; j < num->GetYaxis()->GetNbins(); ++j) {
      for (int k = 0; k < num->GetZaxis()->GetNbins(); ++k) {
        double num_frac_error = num->GetBinError(i + 1, j + 1, k + 1) /
                                num->GetBinContent(i + 1, j + 1, k + 1);
        double denom_frac_error = denom->GetBinError(i + 1, j + 1, k + 1) /
                                  denom->GetBinContent(i + 1, j + 1, k + 1);

        if (!std::isnormal(num_frac_error) ||
            !std::isnormal(denom_frac_error)) {
          ratio->SetBinContent(i + 1, j + 1, k + 1, 0);
        } else if (num_frac_error > frac_error_threshold) {
          ratio->SetBinContent(i + 1, j + 1, k + 1, 0);
        } else if (denom_frac_error > frac_error_threshold) {
          ratio->SetBinContent(i + 1, j + 1, k + 1, 1);
        }
      }
    }
  }
}

void ScrubLowStatsBins(TH2D *num, TH2D *denom, TH2D *ratio,
                       double frac_error_threshold) {
  for (int i = 0; i < num->GetXaxis()->GetNbins(); ++i) {
    for (int j = 0; j < num->GetYaxis()->GetNbins(); ++j) {
      double num_frac_error =
          num->GetBinError(i + 1, j + 1) / num->GetBinContent(i + 1, j + 1);
      double denom_frac_error =
          denom->GetBinError(i + 1, j + 1) / denom->GetBinContent(i + 1, j + 1);

      if (!std::isnormal(num_frac_error) || !std::isnormal(denom_frac_error)) {
        ratio->SetBinContent(i + 1, j + 1, 0);
      } else if (num_frac_error > frac_error_threshold) {
        ratio->SetBinContent(i + 1, j + 1, 0);
      } else if (denom_frac_error > frac_error_threshold) {
        ratio->SetBinContent(i + 1, j + 1, 1);
      }
    }
  }
}

TH2D *Project3DRatio(TH3D *num, TH3D *denom, const char *projconfig,
                     std::string const &name) {
  TH2D *num_proj = dynamic_cast<TH2D *>(num->Project3D(projconfig));
  num_proj->SetDirectory(nullptr);
  TH2D *denom_proj = dynamic_cast<TH2D *>(denom->Project3D(projconfig));
  denom_proj->SetDirectory(nullptr);

  TH2D *rat = dynamic_cast<TH2D *>(num_proj->Clone(name.c_str()));
  rat->Divide(denom_proj);
  rat->SetDirectory(nullptr);

  ScrubLowStatsBins(num_proj, denom_proj, rat, 1.0 / sqrt(50.0));

  return rat;
}

int fakedatarwgen(std::string const &ifile, std::string const &ofile) {

  TFile fin(ifile.c_str());
  if (fin.IsZombie()) {
    std::cout << "Failed to read " << ifile.c_str() << std::endl;
    return 2;
  }

  TFile fout(ofile.c_str(), "RECREATE");
  if (fin.IsZombie()) {
    std::cout << "Failed to write " << ofile.c_str() << std::endl;
    return 2;
  }

  for (std::string const &species : {"numu", "numub", "nue", "nueb"}) {
    for (std::string const &mode :
         {"CCInc", "CC0Pi", "CC1CPi", "CC1Pi0", "CCOther", "NCInc", "NC0Pi",
          "NC1CPi", "NC1Pi0", "NCOther"}) {

      TH3D *neut_nd280_C = dynamic_cast<TH3D *>(
          GetTH1(&fin, "NEUT/ND280/C/" + species + "/EnuPLepThetaLep_" + mode));
      if (neut_nd280_C) { // ND280
        neut_nd280_C->SetDirectory(nullptr);

        TH3D *genie_nd280 = dynamic_cast<TH3D *>(GetTH1(
            &fin, "GENIE/ND280/C/" + species + "/EnuPLepThetaLep_" + mode));
        genie_nd280->SetDirectory(nullptr);

        TH2D *genie_nd280_yx =
            Project3DRatio(genie_nd280, neut_nd280_C, "yx",
                           std::string("genie_nd280_EnuPLep_") + mode);
        TH2D *genie_nd280_zx =
            Project3DRatio(genie_nd280, neut_nd280_C, "zx",
                           std::string("genie_nd280_EnuThetaLep_") + mode);
        TH2D *genie_nd280_yz =
            Project3DRatio(genie_nd280, neut_nd280_C, "yz",
                           std::string("genie_nd280_PLepThetaLep_") + mode);

        fout.WriteTObject(genie_nd280_yx, genie_nd280_yx->GetName());
        fout.WriteTObject(genie_nd280_yz, genie_nd280_yz->GetName());
        fout.WriteTObject(genie_nd280_zx, genie_nd280_zx->GetName());

        TH3D *rat = dynamic_cast<TH3D *>(genie_nd280->Clone(
            (std::string("t2knd_to_nova_C_") + species + "_" + mode).c_str()));
        rat->Divide(neut_nd280_C);
        ScrubLowStatsBins(genie_nd280, neut_nd280_C, rat, 1.0 / sqrt(25));
        rat->SetDirectory(nullptr);
        fout.WriteTObject(rat, rat->GetName());
      }

      TH3D *neut_nd280_O = dynamic_cast<TH3D *>(
          GetTH1(&fin, "NEUT/ND280/O/" + species + "/EnuPLepThetaLep_" + mode));
      if (neut_nd280_O) { // ND280
        neut_nd280_O->SetDirectory(nullptr);

        TH3D *genie_nd280 = dynamic_cast<TH3D *>(GetTH1(
            &fin, "GENIE/ND280/O/" + species + "/EnuPLepThetaLep_" + mode));
        genie_nd280->SetDirectory(nullptr);

        TH2D *genie_nd280_yx =
            Project3DRatio(genie_nd280, neut_nd280_O, "yx",
                           std::string("genie_nd280_EnuPLep_") + mode);
        TH2D *genie_nd280_zx =
            Project3DRatio(genie_nd280, neut_nd280_O, "zx",
                           std::string("genie_nd280_EnuThetaLep_") + mode);
        TH2D *genie_nd280_yz =
            Project3DRatio(genie_nd280, neut_nd280_O, "yz",
                           std::string("genie_nd280_PLepThetaLep_") + mode);

        fout.WriteTObject(genie_nd280_yx, genie_nd280_yx->GetName());
        fout.WriteTObject(genie_nd280_yz, genie_nd280_yz->GetName());
        fout.WriteTObject(genie_nd280_zx, genie_nd280_zx->GetName());

        TH3D *rat = dynamic_cast<TH3D *>(genie_nd280->Clone(
            (std::string("t2knd_to_nova_O_") + species + "_" + mode).c_str()));
        rat->Divide(neut_nd280_O);
        ScrubLowStatsBins(genie_nd280, neut_nd280_O, rat, 1.0 / sqrt(25));
        rat->SetDirectory(nullptr);
        fout.WriteTObject(rat, rat->GetName());
      }

      TH3D *neut_nd280_H = dynamic_cast<TH3D *>(
          GetTH1(&fin, "NEUT/ND280/H/" + species + "/EnuPLepThetaLep_" + mode));
      if (neut_nd280_H) { // ND280
        neut_nd280_H->SetDirectory(nullptr);

        TH3D *genie_nd280 = dynamic_cast<TH3D *>(GetTH1(
            &fin, "GENIE/ND280/H/" + species + "/EnuPLepThetaLep_" + mode));
        genie_nd280->SetDirectory(nullptr);

        TH2D *genie_nd280_yx =
            Project3DRatio(genie_nd280, neut_nd280_H, "yx",
                           std::string("genie_nd280_EnuPLep_") + mode);
        TH2D *genie_nd280_zx =
            Project3DRatio(genie_nd280, neut_nd280_H, "zx",
                           std::string("genie_nd280_EnuThetaLep_") + mode);
        TH2D *genie_nd280_yz =
            Project3DRatio(genie_nd280, neut_nd280_H, "yz",
                           std::string("genie_nd280_PLepThetaLep_") + mode);

        fout.WriteTObject(genie_nd280_yx, genie_nd280_yx->GetName());
        fout.WriteTObject(genie_nd280_yz, genie_nd280_yz->GetName());
        fout.WriteTObject(genie_nd280_zx, genie_nd280_zx->GetName());

        TH3D *rat = dynamic_cast<TH3D *>(genie_nd280->Clone(
            (std::string("t2knd_to_nova_H_") + species + "_" + mode).c_str()));
        rat->Divide(neut_nd280_H);
        ScrubLowStatsBins(genie_nd280, neut_nd280_H, rat, 1.0 / sqrt(25));
        rat->SetDirectory(nullptr);
        fout.WriteTObject(rat, rat->GetName());
      }

      TH3D *neut_novand_plep_C = dynamic_cast<TH3D *>(
          GetTH1(&fin, "NEUT/NOvAND/C/" + species + "/EnuPLepEAvHad_" + mode));
      if (neut_novand_plep_C) { // NOvAND
        neut_novand_plep_C->SetDirectory(nullptr);

        TH3D *genie_novand = dynamic_cast<TH3D *>(GetTH1(
            &fin, "GENIE/NOvAND/C/" + species + "/EnuPLepEAvHad_" + mode));
        genie_novand->SetDirectory(nullptr);

        TH2D *nova_to_t2k_yx =
            Project3DRatio(neut_novand_plep_C, genie_novand, "yx",
                           std::string("genie_novand_EnuPLep_") + mode);
        TH2D *nova_to_t2k_zx =
            Project3DRatio(neut_novand_plep_C, genie_novand, "zx",
                           std::string("genie_novand_EnuEAvHad_") + mode);
        TH2D *nova_to_t2k_yz =
            Project3DRatio(neut_novand_plep_C, genie_novand, "yz",
                           std::string("genie_novand_PLepEAvHad_") + mode);
        fout.WriteTObject(nova_to_t2k_yx, nova_to_t2k_yx->GetName());
        fout.WriteTObject(nova_to_t2k_zx, nova_to_t2k_zx->GetName());
        fout.WriteTObject(nova_to_t2k_yz, nova_to_t2k_yz->GetName());

        TH3D *rat = dynamic_cast<TH3D *>(neut_novand_plep_C->Clone(
            (std::string("nova_to_t2k_plep_C_") + species + "_" + mode)
                .c_str()));
        rat->Divide(genie_novand);
        ScrubLowStatsBins(neut_novand_plep_C, genie_novand, rat,
                          1.0 / sqrt(25));
        rat->SetDirectory(nullptr);
        fout.WriteTObject(rat, rat->GetName());
      }

      TH3D *neut_novand_plep_H = dynamic_cast<TH3D *>(
          GetTH1(&fin, "NEUT/NOvAND/H/" + species + "/EnuPLepEAvHad_" + mode));
      if (neut_novand_plep_H) { // NOvAND
        neut_novand_plep_H->SetDirectory(nullptr);

        TH3D *genie_novand = dynamic_cast<TH3D *>(GetTH1(
            &fin, "GENIE/NOvAND/H/" + species + "/EnuPLepEAvHad_" + mode));
        genie_novand->SetDirectory(nullptr);

        TH2D *nova_to_t2k_yx =
            Project3DRatio(neut_novand_plep_H, genie_novand, "yx",
                           std::string("genie_novand_EnuPLep_") + mode);
        TH2D *nova_to_t2k_zx =
            Project3DRatio(neut_novand_plep_H, genie_novand, "zx",
                           std::string("genie_novand_EnuEAvHad_") + mode);
        TH2D *nova_to_t2k_yz =
            Project3DRatio(neut_novand_plep_H, genie_novand, "yz",
                           std::string("genie_novand_PLepEAvHad_") + mode);
        fout.WriteTObject(nova_to_t2k_yx, nova_to_t2k_yx->GetName());
        fout.WriteTObject(nova_to_t2k_zx, nova_to_t2k_zx->GetName());
        fout.WriteTObject(nova_to_t2k_yz, nova_to_t2k_yz->GetName());

        TH3D *rat = dynamic_cast<TH3D *>(neut_novand_plep_H->Clone(
            (std::string("nova_to_t2k_plep_H_") + species + "_" + mode)
                .c_str()));
        rat->Divide(genie_novand);
        ScrubLowStatsBins(neut_novand_plep_H, genie_novand, rat,
                          1.0 / sqrt(25));
        rat->SetDirectory(nullptr);
        fout.WriteTObject(rat, rat->GetName());
      }

      TH3D *neut_novand_q2_C = dynamic_cast<TH3D *>(
          GetTH1(&fin, "NEUT/NOvAND/C/" + species + "/EnuQ2EAvHad_" + mode));
      if (neut_novand_q2_C) { // NOvAND
        neut_novand_q2_C->SetDirectory(nullptr);

        TH3D *genie_novand = dynamic_cast<TH3D *>(
            GetTH1(&fin, "GENIE/NOvAND/C/" + species + "/EnuQ2EAvHad_" + mode));
        genie_novand->SetDirectory(nullptr);

        TH2D *nova_to_t2k_yx =
            Project3DRatio(neut_novand_q2_C, genie_novand, "yx",
                           std::string("genie_novand_EnuQ2_") + mode);
        TH2D *nova_to_t2k_yz =
            Project3DRatio(neut_novand_q2_C, genie_novand, "yz",
                           std::string("genie_novand_Q2EAvHad_") + mode);
        fout.WriteTObject(nova_to_t2k_yx, nova_to_t2k_yx->GetName());
        fout.WriteTObject(nova_to_t2k_yz, nova_to_t2k_yz->GetName());

        TH3D *rat = dynamic_cast<TH3D *>(neut_novand_q2_C->Clone(
            (std::string("nova_to_t2k_Q2_C_") + species + "_" + mode).c_str()));
        rat->Divide(genie_novand);
        ScrubLowStatsBins(neut_novand_q2_C, genie_novand, rat, 1.0 / sqrt(25));
        rat->SetDirectory(nullptr);
        fout.WriteTObject(rat, rat->GetName());
      }

      TH3D *neut_novand_q2_H = dynamic_cast<TH3D *>(
          GetTH1(&fin, "NEUT/NOvAND/H/" + species + "/EnuQ2EAvHad_" + mode));
      if (neut_novand_q2_H) { // NOvAND
        neut_novand_q2_H->SetDirectory(nullptr);

        TH3D *genie_novand = dynamic_cast<TH3D *>(
            GetTH1(&fin, "GENIE/NOvAND/H/" + species + "/EnuQ2EAvHad_" + mode));
        genie_novand->SetDirectory(nullptr);

        TH2D *nova_to_t2k_yx =
            Project3DRatio(neut_novand_q2_H, genie_novand, "yx",
                           std::string("genie_novand_EnuQ2_") + mode);
        TH2D *nova_to_t2k_yz =
            Project3DRatio(neut_novand_q2_H, genie_novand, "yz",
                           std::string("genie_novand_Q2EAvHad_") + mode);
        fout.WriteTObject(nova_to_t2k_yx, nova_to_t2k_yx->GetName());
        fout.WriteTObject(nova_to_t2k_yz, nova_to_t2k_yz->GetName());

        TH3D *rat = dynamic_cast<TH3D *>(neut_novand_q2_H->Clone(
            (std::string("nova_to_t2k_Q2_H_") + species + "_" + mode).c_str()));
        rat->Divide(genie_novand);
        ScrubLowStatsBins(neut_novand_q2_H, genie_novand, rat, 1.0 / sqrt(25));
        rat->SetDirectory(nullptr);
        fout.WriteTObject(rat, rat->GetName());
      }

      TH3D *neut_novand_ptlep_C = dynamic_cast<TH3D *>(
          GetTH1(&fin, "NEUT/NOvAND/C/" + species + "/EnuPtLepEAvHad_" + mode));
      if (neut_novand_ptlep_C) { // NOvAND
        neut_novand_ptlep_C->SetDirectory(nullptr);

        TH3D *genie_novand = dynamic_cast<TH3D *>(GetTH1(
            &fin, "GENIE/NOvAND/C/" + species + "/EnuPtLepEAvHad_" + mode));
        genie_novand->SetDirectory(nullptr);

        TH2D *nova_to_t2k_yx =
            Project3DRatio(neut_novand_ptlep_C, genie_novand, "yx",
                           std::string("genie_novand_EnuPtLep_") + mode);
        TH2D *nova_to_t2k_yz =
            Project3DRatio(neut_novand_ptlep_C, genie_novand, "yz",
                           std::string("genie_novand_PtLepEAvHad_") + mode);
        fout.WriteTObject(nova_to_t2k_yx, nova_to_t2k_yx->GetName());
        fout.WriteTObject(nova_to_t2k_yz, nova_to_t2k_yz->GetName());

        TH3D *rat = dynamic_cast<TH3D *>(neut_novand_ptlep_C->Clone(
            (std::string("nova_to_t2k_ptlep_C_") + species + "_" + mode)
                .c_str()));
        rat->Divide(genie_novand);
        ScrubLowStatsBins(neut_novand_ptlep_C, genie_novand, rat,
                          1.0 / sqrt(25));
        rat->SetDirectory(nullptr);
        fout.WriteTObject(rat, rat->GetName());
      }

      TH3D *neut_novand_ptlep_H = dynamic_cast<TH3D *>(
          GetTH1(&fin, "NEUT/NOvAND/H/" + species + "/EnuPtLepEAvHad_" + mode));
      if (neut_novand_ptlep_H) { // NOvAND
        neut_novand_ptlep_H->SetDirectory(nullptr);

        TH3D *genie_novand = dynamic_cast<TH3D *>(GetTH1(
            &fin, "GENIE/NOvAND/H/" + species + "/EnuPtLepEAvHad_" + mode));
        genie_novand->SetDirectory(nullptr);

        TH2D *nova_to_t2k_yx =
            Project3DRatio(neut_novand_ptlep_H, genie_novand, "yx",
                           std::string("genie_novand_EnuPtLep_") + mode);
        TH2D *nova_to_t2k_yz =
            Project3DRatio(neut_novand_ptlep_H, genie_novand, "yz",
                           std::string("genie_novand_PtLepEAvHad_") + mode);
        fout.WriteTObject(nova_to_t2k_yx, nova_to_t2k_yx->GetName());
        fout.WriteTObject(nova_to_t2k_yz, nova_to_t2k_yz->GetName());

        TH3D *rat = dynamic_cast<TH3D *>(neut_novand_ptlep_H->Clone(
            (std::string("nova_to_t2k_ptlep_H_") + species + "_" + mode)
                .c_str()));
        rat->Divide(genie_novand);
        ScrubLowStatsBins(neut_novand_ptlep_H, genie_novand, rat,
                          1.0 / sqrt(25));
        rat->SetDirectory(nullptr);
        fout.WriteTObject(rat, rat->GetName());
      }
    }
  }

  fout.Write();
  fout.Close();
  return 0;
}

int main(int argc, char const * argv[]){
  return fakedatarwgen(argv[1], argv[2]);
}
