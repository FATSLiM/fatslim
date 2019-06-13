#!/usr/bin/env bash
SCRIPTPATH="$( cd "$(dirname "$0")" ; pwd -P )"

for i in "$@"
do
case $i in
    --coverage)
    export CYTHON_COVERAGE=1
    shift # past argument=value
    ;;
    *)
    echo "USAGE: $0 [--coverage]"
    exit 0
    ;;
esac
done


cd ${SCRIPTPATH}/../package
python setup.py build_ext -i

cd ${SCRIPTPATH}/../testsuite
if [[ $CYTHON_COVERAGE == 1 ]]; then \
    pytest --cov-report=xml --cov-report=html --cov=fatslim fatslimtest; else \
    pytest fatslimtest -xv
fi

unset CYTHON_COVERAGE