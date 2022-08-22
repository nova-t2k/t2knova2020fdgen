#!/bin/bash

mkdir -p ValidPlots_T2K
cd ValidPlots_T2K

../bin/ValidPlots.exe -i ../FDSValid/FakeDataValid.root \
	--From Tuned_BANFFPre --To Tuned
cd ..

mkdir -p ValidPlots_T2K_NonQE
cd ValidPlots_T2K_NonQE

../bin/ValidPlots.exe -i ../FDSValid/FakeDataValid.root \
	--From Tuned_BANFFPre --To NonQE
cd ..

mkdir -p ValidPlots_Mnv1Pi
cd ValidPlots_Mnv1Pi

../bin/ValidPlots.exe -i ../FDSValid/FakeDataValid.root \
	--From Tuned_BANFFPre --To Mnv1Pi
