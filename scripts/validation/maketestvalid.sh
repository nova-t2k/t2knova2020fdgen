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
                    -M \
                    -o FDSValid/FakeDataValid_C_numu_numub_T2KND_to_NOvA.root \
                    -d ND280/T2KND_to_NOvA/C/numu &

bin/fakedatavalid.exe -i flattrees/t2knova.flattree.NEUT.ND280.CH.numub.root \
                    -F FDSInputs/FakeDataInputs.root \
                    -H config/FakeDataValidConfig.toml \
                    -W T2KND_to_NOvA \
                    -a C \
                    -M \
                    -o FDSValid/FakeDataValid_C_numu_numub_T2KND_to_NOvA.root \
                    -d ND280/T2KND_to_NOvA/C/numub &
wait

bin/fakedatavalid.exe -i flattrees/t2knova.flattree.NEUT.ND280.CH.numu.root \
                    -F FDSInputs/FakeDataInputs.root \
                    -H config/FakeDataValidConfig.toml \
                    -a C \
                    -M \
                    -o FDSValid/FakeDataValid_C_numu_numub_NEUT.root \
                    -d ND280/T2KNDTune/C/numu &

bin/fakedatavalid.exe -i flattrees/t2knova.flattree.GENIE.ND280.CH.numu.root \
                    -F FDSInputs/FakeDataInputs.root \
                    -H config/FakeDataValidConfig.toml \
                    -a C \
                    -M \
                    -o FDSValid/FakeDataValid_C_numu_numub_GENIE.root \
                    -d ND280/NOvATune/C/numu &

wait

bin/fakedatavalid.exe -i flattrees/t2knova.flattree.NEUT.ND280.CH.numub.root \
                    -F FDSInputs/FakeDataInputs.root \
                    -H config/FakeDataValidConfig.toml \
                    -a C \
                    -M \
                    -o FDSValid/FakeDataValid_C_numu_numub_NEUT.root \
                    -d ND280/T2KNDTune/C/numub &

bin/fakedatavalid.exe -i flattrees/t2knova.flattree.GENIE.ND280.CH.numub.root \
                    -F FDSInputs/FakeDataInputs.root \
                    -H config/FakeDataValidConfig.toml \
                    -a C \
                    -M \
                    -o FDSValid/FakeDataValid_C_numu_numub_GENIE.root \
                    -d ND280/NOvATune/C/numub &

wait

hadd -j 4 FDSValid/FakeDataValid_C_numu_numub.root \
			FDSValid/FakeDataValid_C_numu_numub_NEUT.root \
			FDSValid/FakeDataValid_C_numu_numub_GENIE.root \
            FDSValid/FakeDataValid_C_numu_numub_T2KND_to_NOvA.root
