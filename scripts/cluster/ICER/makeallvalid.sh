#!/bin/bash

g++ -I ../../../ fakedatavalid.C -O3 $(root-config --cflags --glibs) -o fakedatavalid.exe
if [[ ! "$?" == "0" ]]; then
        exit;
fi

declare -A TGTEL
TGTEL["NOvAND_CH"]="C H"
TGTEL["ND280_CH"]="C"
TGTEL["ND280_H2O"]="H O"

GENERATORS=( GENIE NEUT )
SPECIES=( numu numub nue nueb )
DETECTORS=( NOvAND ND280 )

declare -A DET_MATS
DET_MATS["NOvAND"]="CH"
DET_MATS["ND280"]="H2O CH"

declare -A MAT_ELEMENTS
MAT_ELEMENTS["CH"]= "C H"
MAT_ELEMENTS["H20"]= "H O"

rm -f FakeDataValid.root
for TGT in ${DET_MATS["ND280"]}; do
  for SPC in ${SPECIES[@]}; do
    ./fakedatavalid.exe t2knova.flattree.NEUT.ND280.${TGT}.${SPC}.root \
        T2KND_To_NOvA \
        FakeDataValid.root ND280/T2KNDTune_To_NOvATune/${TGT}/${SPC}
 
    ./fakedatavalid.exe t2knova.flattree.NEUT.ND280.${TGT}.${SPC}.root \
        T2KND_To_NOvA_Enu \
        FakeDataValid.root ND280/T2KNDTune_To_NOvATune_Enu/${TGT}/${SPC}
    ./fakedatavalid.exe t2knova.flattree.NEUT.ND280.${TGT}.${SPC}.root \
        T2KND_To_NOvA_Q2 \
        FakeDataValid.root ND280/T2KNDTune_To_NOvATune_Q2/${TGT}/${SPC}

    ./fakedatavalid.exe t2knova.flattree.NEUT.ND280.${TGT}.${SPC}.root \
        UnWeighted \
        FakeDataValid.root ND280/T2KNDTune/${TGT}/${SPC}
    ./fakedatavalid.exe t2knova.flattree.GENIE.ND280.${TGT}.${SPC}.root \
        UnWeighted \
        FakeDataValid.root ND280/NOvATune/${TGT}/${SPC}
  
    for ELE in ${MAT_ELEMENTS["${TGT}"]}; do

      ./fakedatavalid.exe t2knova.flattree.NEUT.ND280.${TGT}.${SPC}.root \
         T2KND_To_NOvA \
         FakeDataValid.root ND280/T2KNDTune_To_NOvATune/${TGT}/${ELE}/${SPC} ${ELE}

      ./fakedatavalid.exe t2knova.flattree.NEUT.ND280.${TGT}.${SPC}.root \
         T2KND_To_NOvA_Enu \
         FakeDataValid.root ND280/T2KNDTune_To_NOvATune_Enu/${TGT}/${ELE}/${SPC} ${ELE}

      ./fakedatavalid.exe t2knova.flattree.NEUT.ND280.${TGT}.${SPC}.root \
         T2KND_To_NOvA_Q2 \
         FakeDataValid.root ND280/T2KNDTune_To_NOvATune_Q2/${TGT}/${ELE}/${SPC} ${ELE}

      ./fakedatavalid.exe t2knova.flattree.NEUT.ND280.${TGT}.${SPC}.root \
         UnWeighted \
         FakeDataValid.root ND280/T2KNDTune/${TGT}/${ELE}/${SPC} ${ELE}

      ./fakedatavalid.exe t2knova.flattree.GENIE.ND280.${TGT}.${SPC}.root \
         UnWeighted \
         FakeDataValid.root ND280/NOvATune/${TGT}/${ELE}/${SPC} ${ELE}

    done

  done
done
 
