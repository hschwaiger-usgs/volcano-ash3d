#!/bin/bash
pushd ../tests/test_01
bash clean.sh
bash run_test01.sh
popd

pushd ../tests/test_02
bash clean.sh
bash run_test02.sh
popd

pushd ../tests/test_03
bash clean.sh
bash run_test03.sh
popd

pushd ../tests/test_04
bash clean.sh
bash run_test04.sh
popd



