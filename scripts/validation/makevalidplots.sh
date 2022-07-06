#!/bin/bash

FDSDENOMTYPE=Tuned
# FDSDENOMTYPE=Generated

mkdir -p ValidPlots_${FDSDENOMTYPE}
cd ValidPlots_${FDSDENOMTYPE}

../bin/ValidPlots.exe -i ../FDSValid_From${FDSDENOMTYPE}/FakeDataValid.root \
	--From ${FDSDENOMTYPE}