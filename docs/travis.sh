#!/bin/bash

if [ ${TRAVIS_OS_NAME} = "linux" ]; then
    julia --project=. --color=yes -e 'using Pkg; Pkg.instantiate();';
    julia --project=. --color=yes docs/make.jl
fi
