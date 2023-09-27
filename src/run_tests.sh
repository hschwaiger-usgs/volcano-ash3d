#!/bin/bash
pushd ../tests/test_01
./clean.sh
./run_test01.sh
popd

pushd ../tests/test_02
./clean.sh
./run_test02.sh
popd

pushd ../tests/test_03
./clean.sh
./run_test03.sh
popd

pushd ../tests/test_04
./clean.sh
./run_test04.sh
popd



