#!/bin/bash

FDSDENOMTYPE=Tuned_BANFFPre
# FDSDENOMTYPE=Tuned_BANFFPost
# FDSDENOMTYPE=Generated

mkdir -p ValidPlots_${FDSDENOMTYPE}
cd ValidPlots_${FDSDENOMTYPE}

../bin/ValidPlots.exe -i ../FDSValid_From${FDSDENOMTYPE}/FakeDataValid.root \
	--From ${FDSDENOMTYPE} --To Tuned
cd ..

mkdir -p ValidPlots_${FDSDENOMTYPE}_NonQE
cd ValidPlots_${FDSDENOMTYPE}_NonQE

../bin/ValidPlots.exe -i ../FDSValid_From${FDSDENOMTYPE}/FakeDataValid.root \
	--From ${FDSDENOMTYPE} --To NonQE
cd ..

# mkdir -p ValidPlots_${FDSDENOMTYPE}_Mnv1Pi
# cd ValidPlots_${FDSDENOMTYPE}_Mnv1Pi

# ../bin/ValidPlots.exe -i ../FDSValid_From${FDSDENOMTYPE}/FakeDataValid.root \
# 	--From ${FDSDENOMTYPE} --To Mnv1Pi