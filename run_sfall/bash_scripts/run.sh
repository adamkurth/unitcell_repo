#!/usr/bin/env bash

# must be run from CCP4 Applications directory
# /Applications/ccp4-8.0/ 
# ./start 
# cd /Users/adamkurth/Documents/vscode/CXFEL_Image_Analysis/CXFEL/unitcell_project/run_sfall/bash_scripts/
# ./run.sh
#!/usr/bin/env bash

# Check if spacegroup name is provided
if [ -z "$1" ]; then
    echo "No spacegroup name provided"
    exit 1
fi

SPACEGROUP="$1"

# Define directories
PDB_DIR="/Users/adamkurth/Documents/vscode/CXFEL_Image_Analysis/CXFEL/unitcell_project/run_sfall/pdb/pdb_$SPACEGROUP"
TXT_OUTDIR="/Users/adamkurth/Documents/vscode/CXFEL_Image_Analysis/CXFEL/unitcell_project/run_sfall/out/out_$SPACEGROUP"
MTZ_OUTDIR="/Users/adamkurth/Documents/vscode/CXFEL_Image_Analysis/CXFEL/unitcell_project/run_sfall/mtz/mtz_$SPACEGROUP"
HKL_OUTDIR="/Users/adamkurth/Documents/vscode/CXFEL_Image_Analysis/CXFEL/unitcell_project/run_sfall/data/data_$SPACEGROUP"

# Ensure the PDB directory exists
if [ ! -d "$PDB_DIR" ]; then
    echo "Error: PDB directory $PDB_DIR does not exist."
    exit 1
fi

# Create output directories if they do not exist
mkdir -p "$TXT_OUTDIR"
mkdir -p "$MTZ_OUTDIR"
mkdir -p "$HKL_OUTDIR"

# Loop through all PDB files in the specified directory
shopt -s nullglob
for PDB_FILE in "$PDB_DIR"/*.pdb; do
    # Check if the file exists
    if [ ! -f "$PDB_FILE" ]; then
        echo "No PDB files found in $PDB_DIR"
        continue
    fi

    BASENAME=$(basename "$PDB_FILE" .pdb)
    MTZ_FILE="$MTZ_OUTDIR/${BASENAME}.mtz"
    TXT_FILE="$TXT_OUTDIR/${BASENAME}.txt"
    HKL_FILE="$HKL_OUTDIR/${BASENAME}.hkl"

    echo "-----------------------------------"
    echo "Starting to process $PDB_FILE"

    # Run sfall
    echo "Running sfall for $PDB_FILE..."
    gtimeout 2m sfall XYZIN "$PDB_FILE" HKLOUT "$HKL_FILE" <<EOF
MODE SFCALC XYZIN
RESOLUTION 2 20
FORM NGAUSS 5
SYMM $SPACEGROUP
END
EOF
    echo "Finished sfall for $PDB_FILE."

    # Check if MTZ file exists before running mtz2various
    if [ ! -f "$MTZ_FILE" ]; then
        echo "MTZ file $MTZ_FILE does not exist, skipping mtz2various"
        continue
    fi

    # Run mtz2various
    echo "Running mtz2various for $MTZ_FILE..."
    gtimeout 2m mtz2various HKLIN "$MTZ_FILE" HKLOUT "$HKL_FILE" > "$TXT_FILE" <<EOF
LABIN H=H K=K L=L FC=FC PHIC=PHIC
OUTPUT USER '(3I4,2F9.1)'
EOF
    echo "Finished mtz2various for $MTZ_FILE."

    echo "Processed $PDB_FILE"
    echo "-----------------------------------"

done
shopt -u nullglob
