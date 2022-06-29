#include <cmath>
#include <iostream>
#include <string>

#include "ChannelHistCollections.h"
#include "T2KNOvA/FakeDataHelper.hxx"
#include "TFile.h"
#include "TH2.h"
#include "TH3.h"

using namespace t2knova;

bool DoNEUT = true;
bool DoNOvA = true;

bool DoFDS = true;

std::string FromT2KTUNE = "Generated";
std::string ToNOvATUNE = "2020";

std::string ToT2KTUNE = "BANFF_POST";
std::string FromNOvATUNE = "Generated";

std::vector<std::string> FDSToTunes = {"BANFF_PRE", "Mnv1Pi", "NonQE"};

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
    for (t2knova::selection sel_id : AllSelectionList) {
      std::string selection = SelectionList[sel_id];

      if (DoNEUT) {
        for (std::string const &targetnuc : {"C", "H", "O"}) {
          for (std::string const &TOTUNE :
               DoFDS ? std::vector<std::string>{ToNOvATUNE, FDSToTunes[0],
                                                FDSToTunes[1], FDSToTunes[2]}
                     : std::vector<std::string>{ToNOvATUNE}) {

            std::string FromT2KHistName = "NEUT/ND280/" + targetnuc + "/" +
                                          species + "/" + FromT2KTUNE +
                                          "/EnuPLepThetaLep_" + selection;

            std::unique_ptr<TH3> neut_nd280_EnuPLepThetaLep =
                GetTH<TH3>(fin, FromT2KHistName);
            if (neut_nd280_EnuPLepThetaLep) { // ND280
              neut_nd280_EnuPLepThetaLep->SetDirectory(nullptr);

              std::unique_ptr<TH3> genie_nd280 = GetTH<TH3>(
                  fin, "GENIE/ND280/" + targetnuc + "/" + species + "/" +
                           ToNOvATUNE + "/EnuPLepThetaLep_" + selection);
              genie_nd280->SetDirectory(nullptr);

              TDirectory *dir = MakeDirectoryStructure(
                  fout.get(), FromT2KTUNE + "_to_" + ToNOvATUNE +
                                  "/EnuPLepThetaLep/" + targetnuc + "/" +
                                  species + "/");

              std::unique_ptr<TH3> rat(
                  dynamic_cast<TH3 *>(genie_nd280->Clone(selection.c_str())));
              rat->Divide(neut_nd280_EnuPLepThetaLep.get());
              ScrubLowStatsBins(genie_nd280, neut_nd280_EnuPLepThetaLep, rat,
                                MaxFracError);
              rat->SetDirectory(nullptr);
              dir->WriteTObject(rat.get(), rat->GetName());

            } else {
              std::cout << "[WARN]: Expected to find " << FromT2KHistName
                        << std::endl;
            }
          }
        }
      }
    }

    for (t2knova::selection sel_id : AllSelectionList) {
      std::string selection = SelectionList[sel_id];

      if (DoNOvA) {
        for (std::string const &targetnuc : {"C", "H"}) {

          for (std::string const &TOTUNE :
               DoFDS ? std::vector<std::string>{ToT2KTUNE, FDSToTunes[0],
                                                FDSToTunes[1], FDSToTunes[2]}
                     : std::vector<std::string>{ToT2KTUNE}) {

            std::string FromNOvAHistName = "NEUT/NOvAND/" + targetnuc + "/" +
                                           species + "/" + TOTUNE +
                                           "/EnuPtLepEAvHad_" + selection;
            std::unique_ptr<TH3> neut_novand_plep =
                GetTH<TH3>(fin, FromNOvAHistName);
            if (neut_novand_plep) { // NOvAND
              neut_novand_plep->SetDirectory(nullptr);

              std::unique_ptr<TH3> genie_novand = GetTH<TH3>(
                  fin, "GENIE/NOvAND/" + targetnuc + "/" + species + "/" +
                           FromNOvATUNE + "/EnuPtLepEAvHad_" + selection);
              genie_novand->SetDirectory(nullptr);

              TDirectory *dir = MakeDirectoryStructure(
                  fout.get(), FromNOvATUNE + "_to_" + TOTUNE +
                                  "/EnuPtLepEAvHad/" + targetnuc + "/" +
                                  species + "/");

              std::unique_ptr<TH3> rat(dynamic_cast<TH3 *>(
                  neut_novand_plep->Clone(selection.c_str())));
              rat->Divide(genie_novand.get());
              ScrubLowStatsBins(neut_novand_plep, genie_novand, rat,
                                MaxFracError);
              rat->SetDirectory(nullptr);
              dir->WriteTObject(rat.get(), rat->GetName());
            } else {
              std::cout << "[WARN]: Expected to find " << FromNOvAHistName
                        << std::endl;
            }
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
