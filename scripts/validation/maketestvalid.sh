#!/bin/bash

rm -f FDSValid/FakeDataValid_C_numu*.root

mkdir -p FDSValid

set -e
set -x

bin/fakedatavalid.exe -i flattrees/t2knova.flattree.NEUT.ND280.CH.numu.root \
                    -F FDSInputs/FakeDataInputs.root \
                    -H config/FakeDataValidConfig.toml \
                    -W T2KND_to_NOvA \
                    -a C \
                    -o FDSValid/FakeDataValid_C_numu_T2KND_to_NOvA.root \
                    -d ND280/T2KND_to_NOvA/C/numu &

bin/fakedatavalid.exe -i flattrees/t2knova.flattree.NEUT.ND280.CH.numu.root \
                    -F FDSInputs/FakeDataInputs.root \
                    -H config/FakeDataValidConfig.toml \
                    -W T2KND_to_NOvA_EnuKludge \
                    -a C \
                    -o FDSValid/FakeDataValid_C_numu_T2KND_to_NOvA_EnuKludge.root \
                    -d ND280/T2KND_to_NOvA_EnuKludge/C/numu &

bin/fakedatavalid.exe -i flattrees/t2knova.flattree.NEUT.ND280.CH.numu.root \
                    -F FDSInputs/FakeDataInputs.root \
                    -H config/FakeDataValidConfig.toml \
                    -a C \
                    -o FDSValid/FakeDataValid_C_numu_NEUT.root \
                    -d ND280/T2KNDTune/C/numu &

bin/fakedatavalid.exe -i flattrees/t2knova.flattree.GENIE.ND280.CH.numu.root \
                    -F FDSInputs/FakeDataInputs.root \
                    -H config/FakeDataValidConfig.toml \
                    -a C \
                    -o FDSValid/FakeDataValid_C_numu_GENIE.root \
                    -d ND280/NOvATune/C/numu &

wait

hadd -j 4 FDSValid/FakeDataValid_C_numu.root \
			FDSValid/FakeDataValid_C_numu_NEUT.root \
			FDSValid/FakeDataValid_C_numu_GENIE.root \
            FDSValid/FakeDataValid_C_numu_T2KND_to_NOvA.root \
            FDSValid/FakeDataValid_C_numu_T2KND_to_NOvA_EnuKludge.root
