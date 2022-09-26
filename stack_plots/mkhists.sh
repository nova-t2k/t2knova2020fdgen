#!/bin/sh

rm FakeData*.root

# ${det}/${gen}/${LASTDIRNAME}/${tgtel}/${spec}"

../bin/fakedatahists.exe -i ../flattrees/t2knova.flattree.NEUT.NOvAND.CH.BANFF_PRE.numu.root \
  -H ../config/FakeDataConfig.toml \
  -o FakeDataHists_1.root \
  -e NOvAND --FDS Generated -d NOvAND/NEUT/Generated/CH/numu &

  ../bin/fakedatahists.exe -i ../flattrees/t2knova.flattree.NEUT.NOvAND.CH.BANFF_PRE.numu.root \
  -H ../config/FakeDataConfig.toml \
  -o FakeDataHists_2.root \
  -e NOvAND --FDS NDTuned -d NOvAND/NEUT/BANFF_PRE/CH/numu &

../bin/fakedatahists.exe -i ../flattrees/t2knova.flattree.GENIE.NOvAND.CH.2020.numu.root \
  -H ../config/FakeDataConfig.toml \
  -o FakeDataHists_3.root \
  -e NOvAND --FDS Generated -d NOvAND/GENIE/Generated/CH/numu &

  ../bin/fakedatahists.exe -i ../flattrees/t2knova.flattree.GENIE.NOvAND.CH.2020.numu.root \
  -H ../config/FakeDataConfig.toml \
  -o FakeDataHists_4.root \
  -e NOvAND --FDS NDTuned -d NOvAND/GENIE/2020/CH/numu &



../bin/fakedatavalid.exe -i ../flattrees/t2knova.flattree.NEUT.NOvAND.CH.BANFF_PRE.numu.root \
  -H ../config/FakeDataValidConfig_NOvAND.toml \
  -o FakeDataValid_1.root \
  -e NOvAND -T Generated -d NOvAND/NEUT/Generated/CH/numu &

  ../bin/fakedatavalid.exe -i ../flattrees/t2knova.flattree.NEUT.NOvAND.CH.BANFF_PRE.numu.root \
  -H ../config/FakeDataValidConfig_NOvAND.toml \
  -o FakeDataValid_2.root \
  -e NOvAND -T NDTuned -d NOvAND/NEUT/BANFF_PRE/CH/numu &

../bin/fakedatavalid.exe -i ../flattrees/t2knova.flattree.GENIE.NOvAND.CH.2020.numu.root \
  -H ../config/FakeDataValidConfig_NOvAND.toml \
  -o FakeDataValid_3.root \
  -e NOvAND -T Generated -d NOvAND/GENIE/Generated/CH/numu &

  ../bin/fakedatavalid.exe -i ../flattrees/t2knova.flattree.GENIE.NOvAND.CH.2020.numu.root \
  -H ../config/FakeDataValidConfig_NOvAND.toml \
  -o FakeDataValid_4.root \
  -e NOvAND -T NDTuned -d NOvAND/GENIE/2020/CH/numu &

wait

hadd FakeDataValid.root FakeDataValid_*.root &
hadd FakeDataHists.root FakeDataHists_*.root &

wait 

rm FakeData*_*.root