#!/bin/bash

# mkdir -p ValidPlots_T2K
# cd ValidPlots_T2K

  # ../bin/ValidPlots.exe -i ../FDSValid_FromICER/FakeDataValid.root \
  # 	--From Tuned_BANFFPre --To Tuned

  # ../bin/ValidPlots.exe -i ../FDSValid_FromICER/FakeDataValid.root \
  #   --To Tuned_BANFFPost --NOvAND

# cd ..

mkdir -p ValidPlots_T2K_NonQE
cd ValidPlots_T2K_NonQE

  ../bin/ValidPlots.exe -i ../FDSValid_FromICER/FakeDataValid.root \
  --From Generated --To NonQE

  ../bin/ValidPlots.exe -i ../FDSValid_FromICER/FakeDataValid.root \
    --To NonQE --NOvAND
cd ..

mkdir -p ValidPlots_Mnv1Pi
cd ValidPlots_Mnv1Pi

  ../bin/ValidPlots.exe -i ../FDSValid_FromICER/FakeDataValid.root \
  	--From Tuned_BANFFPre --To Mnv1Pi

  ../bin/ValidPlots.exe -i ../FDSValid_FromICER/FakeDataValid.root \
    --To Mnv1Pi --NOvAND
