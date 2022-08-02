#include <cmath>
#include <iostream>
#include <string>

#include "ChannelHistCollections.h"
#include "T2KNOvA/FakeDataHelper.hxx"
#include "TFile.h"
#include "TH1.h"

using namespace t2knova;

bool DoNEUT = true;
bool DoNEUTToFDS = true;
bool DoNOvA = true;

bool DoFDS = true;

constexpr double NMinEvs = 1;
const double MaxFracError = 1.0 / std::sqrt(NMinEvs);

std::string FromT2KTUNE = "Generated";
std::string ToNOvATUNE = "2020";

std::string FromNOvATUNE = "Generated";
std::string ToT2KTUNE = "BANFF_POST";

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
        for (std::string const &proj : {"EnuPLepThetaLep", "Enu"}) {
          for (std::string const &targetnuc : {"C", "H", "O", "CH", "H2O"}) {
            for (std::string const &TOTUNE :
                 DoFDS ? std::vector<std::string>{ToNOvATUNE}
                       : std::vector<std::string>{ToNOvATUNE}) {

              std::string FromT2KHistName = "NEUT/ND280/" + targetnuc + "/" +
                                            species + "/" + FromT2KTUNE + "/" +
                                            proj + "_" + selection;

              std::unique_ptr<TH1> neut_nd280_EnuPLepThetaLep =
                  GetTH<TH1>(fin, FromT2KHistName);
              std::unique_ptr<TH1> genie_nd280 = GetTH<TH1>(
                  fin, "GENIE/ND280/" + targetnuc + "/" + species + "/" +
                           TOTUNE + "/" + proj + "_" + selection);
              if (neut_nd280_EnuPLepThetaLep && genie_nd280) { // ND280
                neut_nd280_EnuPLepThetaLep->SetDirectory(nullptr);
                genie_nd280->SetDirectory(nullptr);

                TDirectory *dir = MakeDirectoryStructure(
                    fout.get(), FromT2KTUNE + "_to_" + TOTUNE + "/" + proj +
                                    "/" + targetnuc + "/" + species + "/");

                std::unique_ptr<TH1> rat(
                    dynamic_cast<TH1 *>(genie_nd280->Clone(selection.c_str())));
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

            for (std::string const &TOTUNE :
                 DoNEUTToFDS
                     ? std::vector<std::string>{FDSToTunes[0], FDSToTunes[1],
                                                FDSToTunes[2]}
                     : std::vector<std::string>{}) {

              if(TOTUNE == FromT2KTUNE){
                continue;
              }

              std::string FromT2KHistName = "NEUT/ND280/" + targetnuc + "/" +
                                            species + "/" + FromT2KTUNE + "/" +
                                            proj + "_" + selection;

              std::unique_ptr<TH1> neut_nd280_EnuPLepThetaLep =
                  GetTH<TH1>(fin, FromT2KHistName);
              std::unique_ptr<TH1> FDS_nd280 = GetTH<TH1>(
                  fin, "NEUT/ND280/" + targetnuc + "/" + species + "/" +
                           TOTUNE + "/" + proj + "_" + selection);
              if (neut_nd280_EnuPLepThetaLep && FDS_nd280) { // ND280
                neut_nd280_EnuPLepThetaLep->SetDirectory(nullptr);
                FDS_nd280->SetDirectory(nullptr);

                TDirectory *dir = MakeDirectoryStructure(
                    fout.get(), FromT2KTUNE + "_to_" + TOTUNE + "/" + proj +
                                    "/" + targetnuc + "/" + species + "/");

                std::unique_ptr<TH1> rat(
                    dynamic_cast<TH1 *>(FDS_nd280->Clone(selection.c_str())));
                rat->Divide(neut_nd280_EnuPLepThetaLep.get());
                ScrubLowStatsBins(FDS_nd280, neut_nd280_EnuPLepThetaLep, rat,
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
    }

    for (t2knova::selection sel_id : AllSelectionList) {
      std::string selection = SelectionList[sel_id];

      if (DoNOvA) {
        for (std::string const &proj : {"EnuPtLepEAvHad", "Enu"}) {
          for (std::string const &targetnuc : {"C", "H"}) {
            for (std::string const &TOTUNE :
                 DoFDS ? std::vector<std::string>{ToT2KTUNE, FDSToTunes[0],
                                                  FDSToTunes[1], FDSToTunes[2]}
                       : std::vector<std::string>{ToT2KTUNE}) {

              std::string ToNOvAHistName = "NEUT/NOvAND/" + targetnuc + "/" +
                                           species + "/" + TOTUNE + "/" + proj +
                                           "_" + selection;
              std::unique_ptr<TH1> neut_novand_plep =
                  GetTH<TH1>(fin, ToNOvAHistName);
              std::unique_ptr<TH1> genie_novand = GetTH<TH1>(
                  fin, "GENIE/NOvAND/" + targetnuc + "/" + species + "/" +
                           FromNOvATUNE + "/" + proj + "_" + selection);
              if (neut_novand_plep && genie_novand) { // NOvAND
                neut_novand_plep->SetDirectory(nullptr);
                genie_novand->SetDirectory(nullptr);

                TDirectory *dir = MakeDirectoryStructure(
                    fout.get(), FromNOvATUNE + "_to_" + TOTUNE + "/" + proj +
                                    "/" + targetnuc + "/" + species + "/");

                std::unique_ptr<TH1> rat(dynamic_cast<TH1 *>(
                    neut_novand_plep->Clone(selection.c_str())));
                rat->Divide(genie_novand.get());
                ScrubLowStatsBins(neut_novand_plep, genie_novand, rat,
                                  MaxFracError);
                rat->SetDirectory(nullptr);
                dir->WriteTObject(rat.get(), rat->GetName());
              } else {
                std::cout << "[WARN]: Expected to find " << ToNOvAHistName
                          << std::endl;
              }
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

void SayUsage(char const *argv[]) {
  std::cout << "[USAGE]: " << argv[0]
            << "<input hist file> <output file> [opts]"
            << "\n"
               "\t--from-ND280 <tune>      : Denominator tune"
               "\t--to-ND280 <tune>        : Numerator tune"
               "\t--from-NOvAND <tune>     : Denominator tune"
               "\t--to-NOvAND <tune>       : Numerator tune"
            << std::endl;
}

void handleOpts(int argc, char const *argv[]) {
  int opt = 1;
  while (opt < argc) {
    if ((std::string(argv[opt]) == "-?") ||
        (std::string(argv[opt]) == "--help")) {
      SayUsage(argv);
      exit(0);
    } else if (std::string(argv[opt]) == "--from-ND280-NEUT") {
      FromT2KTUNE = argv[++opt];
    } else if (std::string(argv[opt]) == "--to-ND280-GENIE") {
      ToNOvATUNE = argv[++opt];
    } else if (std::string(argv[opt]) == "--from-NOvAND-GENIE") {
      FromNOvATUNE = argv[++opt];
    } else if (std::string(argv[opt]) == "--to-NOvAND-NEUT") {
      ToT2KTUNE = argv[++opt];
    }
    opt++;
  }
}

int main(int argc, char const *argv[]) {
  handleOpts(argc, argv);
  return fakedatarwgen(argv[1], argv[2]);
}
