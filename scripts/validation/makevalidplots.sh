#!/bin/bash

# FDSDENOMTYPE=Tuned_BANFFPre
FDSDENOMTYPE=Tuned_BANFFPost
# FDSDENOMTYPE=Generated

# mkdir -p ValidPlots_T2K_${FDSDENOMTYPE}
# cd ValidPlots_T2K_${FDSDENOMTYPE}

# ../bin/ValidPlots.exe -i ../FDSValid_From${FDSDENOMTYPE}/FakeDataValid.root \
# 	--From ${FDSDENOMTYPE} --To Tuned
# cd ..

mkdir -p ValidPlots_T2K_${FDSDENOMTYPE}_NonQE
cd ValidPlots_T2K_${FDSDENOMTYPE}_NonQE

../bin/ValidPlots.exe -i ../FDSValid_From${FDSDENOMTYPE}/FakeDataValid.root \
	--From ${FDSDENOMTYPE} --To NonQE
cd ..

# mkdir -p ValidPlots_NOvAND_${FDSDENOMTYPE}
# cd ValidPlots_NOvAND_${FDSDENOMTYPE}

# ../bin/ValidPlots.exe -i ../FDSValid_FromTuned/FakeDataValid.root \
# 	--From Tuned --NOvAND
# cd ..

# mkdir -p ValidPlots_NOvAND_Tuned_NonQE
# cd ValidPlots_NOvAND_Tuned_NonQE

# ../bin/ValidPlots.exe -i ../FDSValid_FromTuned/FakeDataValid.root \
# 	--From Tuned --NOvAND --To NonQE
# cd ..


# mkdir -p ValidPlots_NOvAND_Tuned_Mnv1Pi
# cd ValidPlots_NOvAND_Tuned_Mnv1Pi

# ../bin/ValidPlots.exe -i ../FDSValid_FromTuned/FakeDataValid.root \
# 	--From Tuned --NOvAND --To Mnv1Pi
# cd ..


# mkdir -p ValidPlots_${FDSDENOMTYPE}_Mnv1Pi
# cd ValidPlots_${FDSDENOMTYPE}_Mnv1Pi

# ../bin/ValidPlots.exe -i ../FDSValid_From${FDSDENOMTYPE}/FakeDataValid.root \
# 	--From ${FDSDENOMTYPE} --To Mnv1Pi
