#!/bin/bash

if [ ${TRAVIS_OS_NAME} = "linux" ]; then
    julia --project --color=yes -e 'using Pkg; Pkg.instantiate(); import GeometricIntegrators; include(joinpath(dirname(pathof(GeometricIntegrators)), "..", "docs", "make.jl"))';
fi
