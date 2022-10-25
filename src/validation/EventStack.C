#define T2KNOVARW_MERGED_CC0PI

#include "T2KNOvATruthTreeReader.h"

using namespace t2knova;

std::map<selection, std::vector<selection>> MappedModes = {
    {kNOvAFSIMode_CC0Pi, {kCC0pi}},
    {kNOvAFSIMode_CC1cPi, {kCC1cpi}},
    {kNOvAFSIMode_CC1Pi0, {kCC1pi0}},
    {kNOvAFSIMode_CCMultiPi, {kCCmultipi}},
    {kNOvAFSIMode_CCOth, {kCCOther}},
    {kNOvAFSIMode_NCInc, {kNC0pi, kNC1cpi, kNC1pi0, kNCmultipi, kNCOther}},
};

int main(int argc, char const *argv[]) {
  TFile fin(argv[1]);
  if (fin.IsZombie()) {
    std::cout << "Failed to read " << argv[1] << std::endl;
    return 2;
  }

  TTreeReader ttrdr("T2KNOvATruthTree", &fin);

  T2KNOvATruthTreeReader rdr(ttrdr);
  int ent = 0, count = 0;
  while (ttrdr.Next()) {
    auto sels = rdr.GetSelections();

    bool have_NOvA_Mode = MappedModes.count(rdr.NOvAFSIMode());
    bool has_primary = rdr.GetPrimarySelection() != kNoPrimarySel;
    bool primary_matches_NOvA =
        have_NOvA_Mode && has_primary &&
        (std::find(MappedModes[rdr.NOvAFSIMode()].begin(),
                   MappedModes[rdr.NOvAFSIMode()].end(),
                   rdr.GetPrimarySelection()) !=
         MappedModes[rdr.NOvAFSIMode()].end());

    bool is_gamma = (sels[1] == kCC1Gamma) || (sels[1] == kCCNGamma) ||
                    (sels[1] == kNC1Gamma) || (sels[1] == kNCNGamma);

    bool is_CCInc_RW =
        std::find(sels.begin(), sels.end(), kCCInc_RW) != sels.end();
    bool is_CCInc = std::find(sels.begin(), sels.end(), kCCInc) != sels.end();

    bool select = have_NOvA_Mode && has_primary && !primary_matches_NOvA;
    select = is_CCInc && !is_CCInc_RW;

    bool say = select && (count < 100);

    if (say) {
      std::cout << "Event: " << ent << "(Selected: " << count
                << "), Mode = " << rdr.Mode() << std::endl;
      std::cout << rdr.PrintStack() << std::endl;
      std::cout << "- Selected as: " << std::endl;
      for (auto sel : sels) {
        std::cout << " " << (sel == rdr.GetPrimarySelection() ? "* " : "  ")
                  << sel << ": " << SelectionList[sel] << std::endl;
      }
      std::cout << "- NOvAFSIMode: " << rdr.NOvAFSIMode() << std::endl;

      std::cout << "\n-------------\n\n";
    }
    if (select) {
      count++;
    }
    ent++;
  }

  std::cout << "selected " << count << "/" << ent << std::endl;
}