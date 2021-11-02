#include <cmath>
#include <iostream>
#include <string>

#include "ChannelHistCollections.h"
#include "T2KNOvAFakeDataHelper.hxx"
#include "TFile.h"
#include "TH2.h"
#include "TH3.h"

using namespace t2knova;

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
    for (std::string const &targetnuc : {"C", "H", "O"}) {
      for (std::string const &selection : SelectionList) {
        std::string HistSuffix = targetnuc + "_" + species + "_" + selection;

        TH3 *neut_nd280_EnuPLepThetaLep = dynamic_cast<TH3 *>(
            GetTH1(&fin, "NEUT/ND280/" + targetnuc + "/" + species +
                             "/EnuPLepThetaLep_" + selection));
        if (neut_nd280_EnuPLepThetaLep) {  // ND280
          neut_nd280_EnuPLepThetaLep->SetDirectory(nullptr);

          TH3 *genie_nd280 = dynamic_cast<TH3 *>(
              GetTH1(&fin, "GENIE/ND280/" + targetnuc + "/" + species +
                               "/EnuPLepThetaLep_" + selection));
          genie_nd280->SetDirectory(nullptr);

          TH2 *genie_nd280_yx =
              Project3DRatio(genie_nd280, neut_nd280_EnuPLepThetaLep, "yx",
                             std::string("genie_nd280_EnuPLep_") + HistSuffix,
                             1.0 / sqrt(50.0));
          TH2 *genie_nd280_zx = Project3DRatio(
              genie_nd280, neut_nd280_EnuPLepThetaLep, "zx",
              std::string("genie_nd280_EnuThetaLep_") + HistSuffix,
              1.0 / sqrt(50.0));
          TH2 *genie_nd280_yz = Project3DRatio(
              genie_nd280, neut_nd280_EnuPLepThetaLep, "yz",
              std::string("genie_nd280_PLepThetaLep_") + HistSuffix,
              1.0 / sqrt(50.0));

          fout.WriteTObject(genie_nd280_yx, genie_nd280_yx->GetName());
          fout.WriteTObject(genie_nd280_yz, genie_nd280_yz->GetName());
          fout.WriteTObject(genie_nd280_zx, genie_nd280_zx->GetName());

          TH3 *rat = dynamic_cast<TH3 *>(genie_nd280->Clone(
              (std::string("t2knd_to_nova_") + HistSuffix).c_str()));
          rat->Divide(neut_nd280_EnuPLepThetaLep);
          ScrubLowStatsBins(genie_nd280, neut_nd280_EnuPLepThetaLep, rat,
                            1.0 / sqrt(25));
          rat->SetDirectory(nullptr);
          fout.WriteTObject(rat, rat->GetName());
        }

        TH1 *neut_nd280_Enu = dynamic_cast<TH1 *>(GetTH1(
            &fin,
            "NEUT/ND280/" + targetnuc + "/" + species + "/Enu_" + selection));
        if (neut_nd280_Enu) {  // ND280
          neut_nd280_Enu->SetDirectory(nullptr);

          TH1 *genie_nd280 = dynamic_cast<TH1 *>(
              GetTH1(&fin, "GENIE/ND280/" + targetnuc + "/" + species +
                               "/Enu_" + selection));
          genie_nd280->SetDirectory(nullptr);

          TH1 *rat = dynamic_cast<TH1 *>(genie_nd280->Clone(
              (std::string("t2knd_to_nova_Enu_") + HistSuffix).c_str()));
          rat->Divide(neut_nd280_Enu);
          ScrubLowStatsBins(genie_nd280, neut_nd280_Enu, rat, 1.0 / sqrt(25));
          rat->SetDirectory(nullptr);
          fout.WriteTObject(rat, rat->GetName());
        }

        TH1 *neut_nd280_Q2 =
            dynamic_cast<TH1 *>(GetTH1(&fin, "NEUT/ND280/" + targetnuc + "/" +
                                                 species + "/Q2_" + selection));
        if (neut_nd280_Q2) {  // ND280
          neut_nd280_Q2->SetDirectory(nullptr);

          TH1 *genie_nd280 = dynamic_cast<TH1 *>(GetTH1(
              &fin,
              "GENIE/ND280/" + targetnuc + "/" + species + "/Q2_" + selection));
          genie_nd280->SetDirectory(nullptr);

          TH1 *rat = dynamic_cast<TH1 *>(genie_nd280->Clone(
              (std::string("t2knd_to_nova_Q2_") + HistSuffix).c_str()));
          rat->Divide(neut_nd280_Q2);
          ScrubLowStatsBins(genie_nd280, neut_nd280_Q2, rat, 1.0 / sqrt(25));
          rat->SetDirectory(nullptr);
          fout.WriteTObject(rat, rat->GetName());
        }

        TH3 *neut_novand_plep = dynamic_cast<TH3 *>(
            GetTH1(&fin, "NEUT/NOvAND/" + targetnuc + "/" + species +
                             "/EnuPLepEAvHad_" + selection));
        if (neut_novand_plep) {  // NOvAND
          neut_novand_plep->SetDirectory(nullptr);

          TH3 *genie_novand = dynamic_cast<TH3 *>(
              GetTH1(&fin, "GENIE/NOvAND/" + targetnuc + "/" + species +
                               "/EnuPLepEAvHad_" + selection));
          genie_novand->SetDirectory(nullptr);

          TH2 *nova_to_t2k_yx =
              Project3DRatio(neut_novand_plep, genie_novand, "yx",
                             std::string("genie_novand_EnuPLep_") + HistSuffix,
                             1.0 / sqrt(50.0));
          TH2 *nova_to_t2k_zx = Project3DRatio(
              neut_novand_plep, genie_novand, "zx",
              std::string("genie_novand_EnuEAvHad_") + HistSuffix,
              1.0 / sqrt(50.0));
          TH2 *nova_to_t2k_yz = Project3DRatio(
              neut_novand_plep, genie_novand, "yz",
              std::string("genie_novand_PLepEAvHad_") + HistSuffix,
              1.0 / sqrt(50.0));
          fout.WriteTObject(nova_to_t2k_yx, nova_to_t2k_yx->GetName());
          fout.WriteTObject(nova_to_t2k_zx, nova_to_t2k_zx->GetName());
          fout.WriteTObject(nova_to_t2k_yz, nova_to_t2k_yz->GetName());

          TH3 *rat = dynamic_cast<TH3 *>(neut_novand_plep->Clone(
              (std::string("nova_to_t2k_plep_") + HistSuffix).c_str()));
          rat->Divide(genie_novand);
          ScrubLowStatsBins(neut_novand_plep, genie_novand, rat,
                            1.0 / sqrt(25));
          rat->SetDirectory(nullptr);
          fout.WriteTObject(rat, rat->GetName());
        }

        TH3 *neut_novand_q2 = dynamic_cast<TH3 *>(
            GetTH1(&fin, "NEUT/NOvAND/" + targetnuc + "/" + species +
                             "/EnuQ2EAvHad_" + selection));
        if (neut_novand_q2) {  // NOvAND
          neut_novand_q2->SetDirectory(nullptr);

          TH3 *genie_novand = dynamic_cast<TH3 *>(
              GetTH1(&fin, "GENIE/NOvAND/" + targetnuc + "/" + species +
                               "/EnuQ2EAvHad_" + selection));
          genie_novand->SetDirectory(nullptr);

          TH2 *nova_to_t2k_yx =
              Project3DRatio(neut_novand_q2, genie_novand, "yx",
                             std::string("genie_novand_EnuQ2_") + HistSuffix,
                             1.0 / sqrt(50.0));
          TH2 *nova_to_t2k_yz =
              Project3DRatio(neut_novand_q2, genie_novand, "yz",
                             std::string("genie_novand_Q2EAvHad_") + HistSuffix,
                             1.0 / sqrt(50.0));
          fout.WriteTObject(nova_to_t2k_yx, nova_to_t2k_yx->GetName());
          fout.WriteTObject(nova_to_t2k_yz, nova_to_t2k_yz->GetName());

          TH3 *rat = dynamic_cast<TH3 *>(neut_novand_q2->Clone(
              (std::string("nova_to_t2k_Q2_") + HistSuffix).c_str()));
          rat->Divide(genie_novand);
          ScrubLowStatsBins(neut_novand_q2, genie_novand, rat, 1.0 / sqrt(25));
          rat->SetDirectory(nullptr);
          fout.WriteTObject(rat, rat->GetName());
        }

        TH3 *neut_novand_ptlep = dynamic_cast<TH3 *>(
            GetTH1(&fin, "NEUT/NOvAND/" + targetnuc + "/" + species +
                             "/EnuPtLepEAvHad_" + selection));
        if (neut_novand_ptlep) {  // NOvAND
          neut_novand_ptlep->SetDirectory(nullptr);

          TH3 *genie_novand = dynamic_cast<TH3 *>(
              GetTH1(&fin, "GENIE/NOvAND/" + targetnuc + "/" + species +
                               "/EnuPtLepEAvHad_" + selection));
          genie_novand->SetDirectory(nullptr);

          TH2 *nova_to_t2k_yx =
              Project3DRatio(neut_novand_ptlep, genie_novand, "yx",
                             std::string("genie_novand_EnuPtLep_") + HistSuffix,
                             1.0 / sqrt(50.0));
          TH2 *nova_to_t2k_yz = Project3DRatio(
              neut_novand_ptlep, genie_novand, "yz",
              std::string("genie_novand_PtLepEAvHad_") + HistSuffix,
              1.0 / sqrt(50.0));
          fout.WriteTObject(nova_to_t2k_yx, nova_to_t2k_yx->GetName());
          fout.WriteTObject(nova_to_t2k_yz, nova_to_t2k_yz->GetName());

          TH3 *rat = dynamic_cast<TH3 *>(neut_novand_ptlep->Clone(
              (std::string("nova_to_t2k_ptlep_") + HistSuffix).c_str()));
          rat->Divide(genie_novand);
          ScrubLowStatsBins(neut_novand_ptlep, genie_novand, rat,
                            1.0 / sqrt(25));
          rat->SetDirectory(nullptr);
          fout.WriteTObject(rat, rat->GetName());
        }
      }
    }
  }

  fout.Write();
  fout.Close();
  return 0;
}

int main(int argc, char const *argv[]) {
  return fakedatarwgen(argv[1], argv[2]);
}
