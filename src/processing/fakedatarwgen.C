#include <cmath>
#include <iostream>
#include <string>

#include "ChannelHistCollections.h"
#include "T2KNOvA/FakeDataHelper.hxx"
#include "TFile.h"
#include "TH2.h"
#include "TH3.h"

using namespace t2knova;

bool DoExtra = false;
bool DoNEUT = true;
bool DoNOvA = true;

int fakedatarwgen(std::string const &ifile, std::string const &ofile) {
  std::unique_ptr<TFile> fin(new TFile(ifile.c_str()));
  if (fin->IsZombie()) {
    std::cout << "Failed to read " << ifile.c_str() << std::endl;
    return 2;
  }

  std::unique_ptr<TFile> fout(new TFile(ofile.c_str(), "RECREATE"));
  if (fin->IsZombie()) {
    std::cout << "Failed to write " << ofile.c_str() << std::endl;
    return 2;
  }

  for (std::string const &species : {"numu", "numub", "nue", "nueb"}) {
    for (t2knova::selection sel_id : ReWeightSelectionList) {
      std::string selection = SelectionList[sel_id];

      if (DoNEUT) {

        for (std::string const &targetnuc : {"C", "H", "O"}) {

          std::unique_ptr<TH3> neut_nd280_EnuPLepThetaLep =
              GetTH<TH3>(fin, "NEUT/ND280/" + targetnuc + "/" + species +
                                  "/EnuPLepThetaLep_" + selection);
          if (neut_nd280_EnuPLepThetaLep) { // ND280
            neut_nd280_EnuPLepThetaLep->SetDirectory(nullptr);

            std::unique_ptr<TH3> genie_nd280 =
                GetTH<TH3>(fin, "GENIE/ND280/" + targetnuc + "/" + species +
                                    "/EnuPLepThetaLep_" + selection);
            genie_nd280->SetDirectory(nullptr);

            TDirectory *dir = MakeDirectoryStructure(
                fout.get(), "t2knd_to_nova/EnuPLepThetaLep/" + targetnuc + "/" +
                                species + "/");

            if (DoExtra) {
              std::unique_ptr<TH2> yx = Project3DRatio(
                  genie_nd280, neut_nd280_EnuPLepThetaLep, "yx",
                  std::string("EnuPLep_") + selection, MaxFracError);
              std::unique_ptr<TH2> zx = Project3DRatio(
                  genie_nd280, neut_nd280_EnuPLepThetaLep, "zx",
                  std::string("EnuThetaLep_") + selection, MaxFracError);
              std::unique_ptr<TH2> yz = Project3DRatio(
                  genie_nd280, neut_nd280_EnuPLepThetaLep, "yz",
                  std::string("PLepThetaLep_") + selection, MaxFracError);

              dir->WriteTObject(yx.get(), yx->GetName());
              dir->WriteTObject(yz.get(), yz->GetName());
              dir->WriteTObject(zx.get(), zx->GetName());
            }

            std::unique_ptr<TH3> rat(
                dynamic_cast<TH3 *>(genie_nd280->Clone(selection.c_str())));
            rat->Divide(neut_nd280_EnuPLepThetaLep.get());
            ScrubLowStatsBins(genie_nd280, neut_nd280_EnuPLepThetaLep, rat,
                              MaxFracError);
            rat->SetDirectory(nullptr);
            dir->WriteTObject(rat.get(), rat->GetName());

            // Calculate Enu Correction
            std::unique_ptr<TH3> IntegralCheck3D(
                dynamic_cast<TH3 *>(neut_nd280_EnuPLepThetaLep->Clone(
                    (std::string("IntegralCheck3D_t2knd_to_nova_") + selection)
                        .c_str())));

            IntegralCheck3D->Multiply(rat.get());

            std::vector<double> EnuBinning =
                GetTAxisBinning(IntegralCheck3D->GetXaxis());

            std::unique_ptr<TH1> EnuNEUT(
                new TH1F((std::string("Enu_NEUT_") + selection).c_str(),
                         ";E_{#nu} (GeV);d#sigma/dE_{#nu}",
                         EnuBinning.size() - 1, EnuBinning.data()));
            std::unique_ptr<TH1> EnuGENIE(
                new TH1F((std::string("Enu_GENIE_") + selection).c_str(),
                         ";E_{#nu} (GeV);d#sigma/dE_{#nu}",
                         EnuBinning.size() - 1, EnuBinning.data()));
            std::unique_ptr<TH1> MissingPSENuCorrection(new TH1F(
                (std::string("MissingPSENuCorrection_") + selection).c_str(),
                ";E_{#nu} (GeV);d#sigma/dE_{#nu}", EnuBinning.size() - 1,
                EnuBinning.data()));

            EnuNEUT->SetDirectory(nullptr);
            EnuGENIE->SetDirectory(nullptr);
            MissingPSENuCorrection->SetDirectory(nullptr);

            for (int i = 0;
                 i < neut_nd280_EnuPLepThetaLep->GetXaxis()->GetNbins(); ++i) {

              std::unique_ptr<TH2> PLepThetaLep_NEUT =
                  ProjectYZSlice(neut_nd280_EnuPLepThetaLep, i + 1, i + 2);
              std::unique_ptr<TH2> PLepThetaLep_GENIE =
                  ProjectYZSlice(genie_nd280, i + 1, i + 2);
              std::unique_ptr<TH2> PLepThetaLep_NEUT_to_GENIE =
                  ProjectYZSlice(IntegralCheck3D, i + 1, i + 2);

              double int_NEUT = IntegralTH2(PLepThetaLep_NEUT);
              double int_GENIE = IntegralTH2(PLepThetaLep_GENIE);
              double int_NEUT_TO_GENIE =
                  IntegralTH2(PLepThetaLep_NEUT_to_GENIE);

              // std::cout
              //     << "Enu Bin: "
              //     << neut_nd280_EnuPLepThetaLep->GetXaxis()->GetBinLowEdge(i
              //     + 1)
              //     << " "
              //     << neut_nd280_EnuPLepThetaLep->GetXaxis()->GetBinUpEdge(i +
              //     1)
              //     << std::endl;
              // std::cout << "\tNEUT Integral = " << int_NEUT << std::endl;
              // std::cout << "\tGENIE Integral = " << int_GENIE << std::endl;
              // std::cout << "\tNEUT * (GENIE/NEUT) Integral = "
              //           << int_NEUT_TO_GENIE << std::endl;
              // std::cout << "\tNEUT * (GENIE/NEUT) Correction = "
              //           << int_GENIE / int_NEUT_TO_GENIE << std::endl;

              double ratio = std::isnormal(int_NEUT_TO_GENIE)
                                 ? int_GENIE / int_NEUT_TO_GENIE
                                 : 1;

              EnuNEUT->SetBinContent(i + 1, int_NEUT);
              EnuGENIE->SetBinContent(i + 1, int_GENIE);
              MissingPSENuCorrection->SetBinContent(i + 1, ratio);
            }

            dir->WriteTObject(EnuNEUT.get(), EnuNEUT->GetName());
            dir->WriteTObject(EnuGENIE.get(), EnuGENIE->GetName());
            dir->WriteTObject(MissingPSENuCorrection.get(),
                              MissingPSENuCorrection->GetName());
          } else {
            std::cout << "[WARN]: Expected to find "
                      << "NEUT/ND280/" + targetnuc + "/" + species +
                             "/EnuPLepThetaLep_" + selection
                      << std::endl;
          }

          std::unique_ptr<TH1> neut_nd280_Enu =
              GetTH1(fin, "NEUT/ND280/" + targetnuc + "/" + species + "/Enu_" +
                              selection);
          if (neut_nd280_Enu) { // ND280
            neut_nd280_Enu->SetDirectory(nullptr);

            std::unique_ptr<TH1> genie_nd280 =
                GetTH1(fin, "GENIE/ND280/" + targetnuc + "/" + species +
                                "/Enu_" + selection);
            genie_nd280->SetDirectory(nullptr);

            std::unique_ptr<TH1> rat(
                dynamic_cast<TH1 *>(genie_nd280->Clone(selection.c_str())));
            rat->Divide(neut_nd280_Enu.get());
            ScrubLowStatsBins(genie_nd280, neut_nd280_Enu, rat, MaxFracError);
            rat->SetDirectory(nullptr);

            TDirectory *dir = MakeDirectoryStructure(
                fout.get(),
                "t2knd_to_nova/Enu/" + targetnuc + "/" + species + "/");

            dir->WriteTObject(rat.get(), rat->GetName());
          } else {
            std::cout << "[WARN]: Expected to find "
                      << "NEUT/ND280/" + targetnuc + "/" + species + "/Enu_" +
                             selection
                      << std::endl;
          }

          std::unique_ptr<TH1> neut_nd280_Q2 =
              GetTH1(fin, "NEUT/ND280/" + targetnuc + "/" + species + "/Q2_" +
                              selection);
          if (neut_nd280_Q2) { // ND280
            neut_nd280_Q2->SetDirectory(nullptr);

            std::unique_ptr<TH1> genie_nd280 =
                GetTH1(fin, "GENIE/ND280/" + targetnuc + "/" + species +
                                "/Q2_" + selection);
            genie_nd280->SetDirectory(nullptr);

            std::unique_ptr<TH1> rat(
                dynamic_cast<TH1 *>(genie_nd280->Clone(selection.c_str())));
            rat->Divide(neut_nd280_Q2.get());
            ScrubLowStatsBins(genie_nd280, neut_nd280_Q2, rat, MaxFracError);
            rat->SetDirectory(nullptr);

            TDirectory *dir = MakeDirectoryStructure(
                fout.get(),
                "t2knd_to_nova/Q2/" + targetnuc + "/" + species + "/");
            dir->WriteTObject(rat.get(), rat->GetName());
          } else {
            std::cout << "[WARN]: Expected to find "
                      << "NEUT/ND280/" + targetnuc + "/" + species + "/Q2_" +
                             selection
                      << std::endl;
          }
        }
      }
    }

    for (t2knova::selection sel_id : t2knova::AllSelectionList) {
      std::string selection = SelectionList[sel_id];

      if (DoNOvA) {

        for (std::string const &targetnuc : {"C", "H"}) {

          std::unique_ptr<TH3> neut_novand_plep =
              GetTH<TH3>(fin, "NEUT/NOvAND/" + targetnuc + "/" + species +
                                  "/EnuPLepEAvHad_" + selection);
          if (neut_novand_plep) { // NOvAND
            neut_novand_plep->SetDirectory(nullptr);

            std::unique_ptr<TH3> genie_novand =
                GetTH<TH3>(fin, "GENIE/NOvAND/" + targetnuc + "/" + species +
                                    "/EnuPLepEAvHad_" + selection);
            genie_novand->SetDirectory(nullptr);

            TDirectory *dir = MakeDirectoryStructure(
                fout.get(),
                "nova_to_t2k/EnuPLepEAvHad/" + targetnuc + "/" + species + "/");

            if (DoExtra) {
              std::unique_ptr<TH2> nova_to_t2k_yx = Project3DRatio(
                  neut_novand_plep, genie_novand, "yx",
                  std::string("EnuPLep_") + selection, MaxFracError);
              std::unique_ptr<TH2> nova_to_t2k_zx = Project3DRatio(
                  neut_novand_plep, genie_novand, "zx",
                  std::string("EnuEAvHad_") + selection, MaxFracError);
              std::unique_ptr<TH2> nova_to_t2k_yz = Project3DRatio(
                  neut_novand_plep, genie_novand, "yz",
                  std::string("PLepEAvHad_") + selection, MaxFracError);
              dir->WriteTObject(nova_to_t2k_yx.get(),
                                nova_to_t2k_yx->GetName());
              dir->WriteTObject(nova_to_t2k_zx.get(),
                                nova_to_t2k_zx->GetName());
              dir->WriteTObject(nova_to_t2k_yz.get(),
                                nova_to_t2k_yz->GetName());
            }
            std::unique_ptr<TH3> rat(dynamic_cast<TH3 *>(
                neut_novand_plep->Clone(selection.c_str())));
            rat->Divide(genie_novand.get());
            ScrubLowStatsBins(neut_novand_plep, genie_novand, rat,
                              MaxFracError);
            rat->SetDirectory(nullptr);
            dir->WriteTObject(rat.get(), rat->GetName());
          } else {
            std::cout << "[WARN]: Expected to find "
                      << "NEUT/NOvAND/" + targetnuc + "/" + species +
                             "/EnuPLepEAvHad_" + selection
                      << std::endl;
          }

          std::unique_ptr<TH3> neut_novand_q2 =
              GetTH<TH3>(fin, "NEUT/NOvAND/" + targetnuc + "/" + species +
                                  "/EnuQ2EAvHad_" + selection);
          if (neut_novand_q2) { // NOvAND
            neut_novand_q2->SetDirectory(nullptr);

            std::unique_ptr<TH3> genie_novand =
                GetTH<TH3>(fin, "GENIE/NOvAND/" + targetnuc + "/" + species +
                                    "/EnuQ2EAvHad_" + selection);
            genie_novand->SetDirectory(nullptr);

            TDirectory *dir = MakeDirectoryStructure(
                fout.get(),
                "nova_to_t2k/EnuQ2EAvHad/" + targetnuc + "/" + species + "/");

            if (DoExtra) {
              std::unique_ptr<TH2> nova_to_t2k_yx = Project3DRatio(
                  neut_novand_q2, genie_novand, "yx",
                  std::string("EnuQ2_") + selection, MaxFracError);
              std::unique_ptr<TH2> nova_to_t2k_yz = Project3DRatio(
                  neut_novand_q2, genie_novand, "yz",
                  std::string("Q2EAvHad_") + selection, MaxFracError);
              dir->WriteTObject(nova_to_t2k_yx.get(),
                                nova_to_t2k_yx->GetName());
              dir->WriteTObject(nova_to_t2k_yz.get(),
                                nova_to_t2k_yz->GetName());
            }

            std::unique_ptr<TH3> rat(
                dynamic_cast<TH3 *>(neut_novand_q2->Clone(selection.c_str())));
            rat->Divide(genie_novand.get());
            ScrubLowStatsBins(neut_novand_q2, genie_novand, rat, MaxFracError);
            rat->SetDirectory(nullptr);
            dir->WriteTObject(rat.get(), rat->GetName());
          } else {
            std::cout << "[WARN]: Expected to find "
                      << "NEUT/NOvAND/" + targetnuc + "/" + species +
                             "/EnuQ2EAvHad_" + selection
                      << std::endl;
          }

          std::unique_ptr<TH3> neut_novand_ptlep =
              GetTH<TH3>(fin, "NEUT/NOvAND/" + targetnuc + "/" + species +
                                  "/EnuPtLepEAvHad_" + selection);
          if (neut_novand_ptlep) { // NOvAND
            neut_novand_ptlep->SetDirectory(nullptr);

            std::unique_ptr<TH3> genie_novand =
                GetTH<TH3>(fin, "GENIE/NOvAND/" + targetnuc + "/" + species +
                                    "/EnuPtLepEAvHad_" + selection);
            genie_novand->SetDirectory(nullptr);

            TDirectory *dir = MakeDirectoryStructure(
                fout.get(), "nova_to_t2k/EnuPtLepEAvHad/" + targetnuc + "/" +
                                species + "/");

            if (DoExtra) {
              std::unique_ptr<TH2> nova_to_t2k_yx = Project3DRatio(
                  neut_novand_ptlep, genie_novand, "yx",
                  std::string("EnuPtLep_") + selection, MaxFracError);
              std::unique_ptr<TH2> nova_to_t2k_yz = Project3DRatio(
                  neut_novand_ptlep, genie_novand, "yz",
                  std::string("PtLepEAvHad_") + selection, MaxFracError);
              dir->WriteTObject(nova_to_t2k_yx.get(),
                                nova_to_t2k_yx->GetName());
              dir->WriteTObject(nova_to_t2k_yz.get(),
                                nova_to_t2k_yz->GetName());
            }

            std::unique_ptr<TH3> rat(dynamic_cast<TH3 *>(
                neut_novand_ptlep->Clone(selection.c_str())));
            rat->Divide(genie_novand.get());
            ScrubLowStatsBins(neut_novand_ptlep, genie_novand, rat,
                              MaxFracError);
            rat->SetDirectory(nullptr);
            dir->WriteTObject(rat.get(), rat->GetName());
          } else {
            std::cout << "[WARN]: Expected to find "
                      << "NEUT/NOvAND/" + targetnuc + "/" + species +
                             "/EnuPtLepEAvHad_" + selection
                      << std::endl;
          }
        }
      }
    }
  }

  fout->Write();
  fout->Close();
  return 0;
}

int main(int argc, char const *argv[]) {
  return fakedatarwgen(argv[1], argv[2]);
}
