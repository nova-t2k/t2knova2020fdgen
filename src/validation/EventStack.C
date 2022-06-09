#include "T2KNOvATruthTreeReader.h"

using namespace t2knova;

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
    if ((std::find(sels.begin(), sels.end(), kCC1Gamma) != sels.end())) {
      std::cout << "Event: " << ent << "(Selected: " << count << "), Mode = " << rdr.Mode() << std::endl;
      std::cout << rdr.PrintStack() << std::endl;
      std::cout << "- Selected as: ";
      for (auto sel : sels) {
        std::cout << " " << sel << ": " << SelectionList[sel];
      }
      std::cout << "\n-------------\n\n";
      count++;
    }
    ent++;

    if (count >= 100) {
      break;
    }
  }
}