#!/bin/bash

pushd test_01
./clean.sh
./clean_plots.sh
popd

pushd test_02
./clean.sh
./clean_plots.sh
popd

pushd test_03
./clean.sh
./clean_plots.sh
popd

pushd test_04
./clean.sh
./clean_plots.sh
popd


