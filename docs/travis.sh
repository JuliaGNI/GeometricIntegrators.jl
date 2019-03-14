#!/bin/bash

if [ ${TRAVIS_OS_NAME} = "linux" ]; then
    julia --color=yes -e 'import GeometricIntegrators; include(joinpath(dirname(pathof(GeometricIntegrators)), "..", "docs", "make.jl"))';
fi
