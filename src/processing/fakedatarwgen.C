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
bool DoFDS = true;
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
          if (neut_nd280_EnuPLepThetaLep) {  // ND280
            neut_nd280_EnuPLepThetaLep->SetDirectory(nullptr);

            std::unique_ptr<TH3> genie_nd280 =
                GetTH<TH3>(fin, "GENIE/ND280/" + targetnuc + "/" + species +
                                    "/EnuPLepThetaLep_" + selection);
            genie_nd280->SetDirectory(nullptr);

            TDirectory *dir = MakeDirectoryStructure(
                fout.get(), "t2knd_to_nova/EnuPLepThetaLep/" + targetnuc + "/" +
                                species + "/");

            std::unique_ptr<TH3> rat(
                dynamic_cast<TH3 *>(genie_nd280->Clone(selection.c_str())));
            rat->Divide(neut_nd280_EnuPLepThetaLep.get());
            ScrubLowStatsBins(genie_nd280, neut_nd280_EnuPLepThetaLep, rat,
                              MaxFracError);
            rat->SetDirectory(nullptr);
            dir->WriteTObject(rat.get(), rat->GetName());

          } else {
            std::cout << "[WARN]: Expected to find "
                      << "NEUT/ND280/" + targetnuc + "/" + species +
                             "/EnuPLepThetaLep_" + selection
                      << std::endl;
          }
        }
      }
    }

    for (t2knova::selection sel_id : ReWeightSelectionList) {
      std::string selection = SelectionList[sel_id];

      if (DoNOvA) {
        for (std::string const &targetnuc : {"C", "H"}) {
          std::unique_ptr<TH3> neut_novand_plep =
              GetTH<TH3>(fin, "NEUT/NOvAND/" + targetnuc + "/" + species +
                                  "/EnuPtLepEAvHad_" + selection);
          if (neut_novand_plep) {  // NOvAND
            neut_novand_plep->SetDirectory(nullptr);

            std::unique_ptr<TH3> genie_novand =
                GetTH<TH3>(fin, "GENIE/NOvAND/" + targetnuc + "/" + species +
                                    "/EnuPtLepEAvHad_" + selection);
            genie_novand->SetDirectory(nullptr);

            TDirectory *dir = MakeDirectoryStructure(
                fout.get(),
                "nova_to_t2k/EnuPtLepEAvHad/" + targetnuc + "/" + species + "/");

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
