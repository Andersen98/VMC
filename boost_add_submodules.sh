#!/bin/bash
declare -a StringArray=("format" "serialization" "mpi" "math"
"algorithm" "container" "smart_ptr")

for module in "${StringArray[@]}"; do
	git submodule add https://github.com/boostorg/$module.git ext/boost/$module
done