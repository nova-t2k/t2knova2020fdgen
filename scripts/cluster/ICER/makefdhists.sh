#!/bin/bash

rm FakeDataHists.root
g++ fakedatahists.C $(root-config --cflags --glibs) -o fakedatahists.exe
if [[ ! "$?" == "0" ]]; then
	exit;
fi

for gen in GENIE NEUT; do
  for det in NOvAND ND280; do
    for spec in numu numub nue nueb; do

     TUNE="nova2020"
     if [ ${gen} == "NEUT" ]; then
       TUNE="banffpost"
     fi 

     SEL="NOvA"
     if [ ${det} == "ND280" ]; then
       SEL="T2K"
     fi

     BEAMMODE="FHC"
     if [ ${spec} == "numub" ] || [ ${spec} == "nueb" ]; then
       BEAMMODE="RHC"
     fi


    echo ./fakedatahists.exe ../flat.${TUNE}.${det}.${spec}.CH.${BEAMMODE}.${gen}.root ${SEL} FakeDataHists.root ${gen}/${det}/${spec}
    done
  done
done

echo root -l -b -q "fakedatarwgen.C(\"FakeDataHists.root\",\"FakeDataInputs.root\")"
