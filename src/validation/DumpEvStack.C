#include "T2KNOvATruthTreeReader.h"

int main(int argc, char const *argv[]) {
  TFile fin(argv[1]);
  if (fin.IsZombie()) {
    std::cout << "Failed to read " << argv[1] << std::endl;
    return 2;
  }

  TTreeReader ttrdr("T2KNOvATruthTree", &fin);

  t2knova::T2KNOvATruthTreeReader rdr(ttrdr);
  while (ttrdr.Next()) {
    auto sels = rdr.GetSelections();
    if ((rdr.Mode() == -1) && (std::find(sels.begin(), sels.end(),
                                         t2knova::kCCOther) != sels.end())) {
      std::cout << rdr.PrintStack() << std::endl;

      std::cout << "- Selected as: ";
      for (auto sel : sels) {
        std::cout << " " << sel << ": " << t2knova::SelectionList[sel];
      }
      std::cout << "\n";
    }
  }
}