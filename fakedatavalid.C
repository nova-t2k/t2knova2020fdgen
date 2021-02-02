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

void fakedatarwgen(std::string const &ifile, std::string const &ofile) {

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

      TH3D *neut_nd280 = dynamic_cast<TH3D *>(
          GetTH1(&fin, "NEUT/ND280/" + species + "/EnuPLepThetaLep_" + mode));
      if (neut_nd280) { // ND280
        neut_nd280->SetDirectory(nullptr);

        TH3D *genie_nd280 = dynamic_cast<TH3D *>(GetTH1(
            &fin, "GENIE/ND280/" + species + "/EnuPLepThetaLep_" + mode));
        genie_nd280->SetDirectory(nullptr);

        TH2D *genie_nd280_yx =
            Project3DRatio(genie_nd280, neut_nd280, "yx",
                           std::string("genie_nd280_EnuPLep_") + mode);
        TH2D *genie_nd280_zx =
            Project3DRatio(genie_nd280, neut_nd280, "zx",
                           std::string("genie_nd280_EnuThetaLep_") + mode);
        TH2D *genie_nd280_yz =
            Project3DRatio(genie_nd280, neut_nd280, "yz",
                           std::string("genie_nd280_PLepThetaLep_") + mode);

        fout.WriteTObject(genie_nd280_yx, genie_nd280_yx->GetName());
        fout.WriteTObject(genie_nd280_yz, genie_nd280_yz->GetName());
        fout.WriteTObject(genie_nd280_zx, genie_nd280_zx->GetName());

        TH3D *rat = dynamic_cast<TH3D *>(genie_nd280->Clone(
            (std::string("t2knd_to_nova_") + species + "_" + mode).c_str()));
        rat->Divide(neut_nd280);
        ScrubLowStatsBins(genie_nd280, neut_nd280, rat, 1.0 / sqrt(25));
        rat->SetDirectory(nullptr);
        fout.WriteTObject(rat, rat->GetName());
      }

      TH3D *neut_sk = dynamic_cast<TH3D *>(
          GetTH1(&fin, "NEUT/SK/" + species + "/EnuPLepThetaLep_" + mode));
      if (neut_sk) { // SK
        neut_sk->SetDirectory(nullptr);

        TH3D *genie_sk = dynamic_cast<TH3D *>(
            GetTH1(&fin, "GENIE/SK/" + species + "/EnuPLepThetaLep_" + mode));
        genie_sk->SetDirectory(nullptr);

        TH2D *genie_sk_yx = Project3DRatio(
            genie_sk, neut_sk, "yx", std::string("genie_sk_EnuPLep_") + mode);
        TH2D *genie_sk_zx =
            Project3DRatio(genie_sk, neut_sk, "zx",
                           std::string("genie_sk_EnuThetaLep_") + mode);
        TH2D *genie_sk_yz =
            Project3DRatio(genie_sk, neut_sk, "yz",
                           std::string("genie_sk_PLepThetaLep_") + mode);
        fout.WriteTObject(genie_sk_yx, genie_sk_yx->GetName());
        fout.WriteTObject(genie_sk_yz, genie_sk_yz->GetName());
        fout.WriteTObject(genie_sk_zx, genie_sk_zx->GetName());

        TH3D *rat = dynamic_cast<TH3D *>(genie_sk->Clone(
            (std::string("sk_to_nova_") + species + "_" + mode).c_str()));
        rat->Divide(neut_sk);
        ScrubLowStatsBins(genie_sk, neut_sk, rat, 1.0 / sqrt(25));
        rat->SetDirectory(nullptr);
        fout.WriteTObject(rat, rat->GetName());
      }

      TH3D *neut_novand = dynamic_cast<TH3D *>(
          GetTH1(&fin, "NEUT/NOvAND/" + species + "/EnuPLepEAvHad_" + mode));
      if (neut_novand) { // NOvAND
        neut_novand->SetDirectory(nullptr);

        TH3D *genie_novand = dynamic_cast<TH3D *>(
            GetTH1(&fin, "GENIE/NOvAND/" + species + "/EnuPLepEAvHad_" + mode));
        genie_novand->SetDirectory(nullptr);

        TH2D *nova_to_t2k_yx =
            Project3DRatio(neut_novand, genie_novand, "yx",
                           std::string("genie_novand_EnuPLep_") + mode);
        TH2D *nova_to_t2k_zx =
            Project3DRatio(neut_novand, genie_novand, "zx",
                           std::string("genie_novand_EnuEAvHad_") + mode);
        TH2D *nova_to_t2k_yz =
            Project3DRatio(neut_novand, genie_novand, "yz",
                           std::string("genie_novand_PLepEAvHad_") + mode);
        fout.WriteTObject(nova_to_t2k_yx, nova_to_t2k_yx->GetName());
        fout.WriteTObject(nova_to_t2k_zx, nova_to_t2k_zx->GetName());
        fout.WriteTObject(nova_to_t2k_yz, nova_to_t2k_yz->GetName());

        TH3D *rat = dynamic_cast<TH3D *>(neut_novand->Clone(
            (std::string("nova_to_t2k_") + species + "_" + mode).c_str()));
        rat->Divide(genie_novand);
        ScrubLowStatsBins(neut_novand, genie_novand, rat, 1.0 / sqrt(25));
        rat->SetDirectory(nullptr);
        fout.WriteTObject(rat, rat->GetName());
      }
    }
  }

  fout.Write();
  fout.Close();
}
