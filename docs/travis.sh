#!/bin/bash

if [ ${TRAVIS_OS_NAME} = "linux" ]; then
    cd docs
    julia --project=. --color=yes -e 'using Pkg; Pkg.instantiate();';
    julia --project=. --color=yes docs/tutorial.jl
    julia --project=. --color=yes docs/make.jl
fi
