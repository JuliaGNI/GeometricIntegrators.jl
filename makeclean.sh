#!/bin/bash

find . -type f -name '*.jl.*.mem' -exec rm {} +
find . -type f -name '.DS_Store' -exec rm {} +
find . -type d -name '.ipynb_checkpoints' -exec rm -Rf {} +
