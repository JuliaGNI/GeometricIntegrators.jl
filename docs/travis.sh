#!/bin/bash

if [ ${TRAVIS_OS_NAME} = "linux" ]; then
    julia --color=yes -e 'Pkg.add("Documenter")';
    julia --color=yes -e 'cd(Pkg.dir("PoincareInvariants")); include(joinpath("docs", "make.jl"))';
fi
