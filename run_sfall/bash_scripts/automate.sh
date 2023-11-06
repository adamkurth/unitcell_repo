#!/usr/bin/env bash

# Check if spacegroup name is provided
if [ -z "$1" ]
then
    echo "No spacegroup name provided"
    exit 1
fi

SPACEGROUP="$1"

# Start CCP4 environment and run your script
/Applications/ccp4-8.0/ ./start /Users/adamkurth/Documents/vscode/CXFEL_Image_Analysis/CXFEL/unitcell_project/run_sfall/bash_scripts/run.sh "$SPACEGROUP"
