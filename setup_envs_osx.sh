#!/bin/bash

# Define the list of environment files
ENV_FILES=(
    "lib/vcftools.yml"
    "lib/align.yml"
    "lib/fpgVcf2FastaEnv.yml"
    "lib/vcfkit.yml"
    "lib/callVar.yml"
    "lib/phylo.yml"
    "lib/fpgDenovo.yml"
    "lib/fpgSNPannotation.yml"
    "lib/mask.yml"
    "lib/fpgtrimReads.yml"
)

echo "=========================================="
echo "Starting Conda Environment Setup for FungiPhyloGen"
echo "=========================================="

# check if conda is installed
if ! command -v conda &> /dev/null; then
    echo "Error: Conda is not installed or not in your PATH."
    exit 1
fi

# Loop through each file
for env_file in "${ENV_FILES[@]}"; do
    if [ -f "$env_file" ]; then
        # OPTIMIZATION 1: Safer extraction of env name (handles spacing and Windows carriage returns)
        env_name=$(awk '/^name:/ {print $2; exit}' "$env_file" | tr -d '\r')
        
        echo "Processing $env_file (Env Name: $env_name)..."

        # OPTIMIZATION 2: Exact word match for environment checking to prevent partial string matches
        if conda info --envs | awk '{print $1}' | grep -Fxq "$env_name"; then
            echo "  [SKIP] Environment '$env_name' already exists."
        else
            echo "  [CREATE] Creating environment '$env_name' using osx-64 emulation..."
            
            # OPTIMIZATION 3: Force Intel architecture for Apple Silicon compatibility
            # Removed the '-y' flag as `conda env create` does not strictly support it and it can cause warnings
            CONDA_SUBDIR=osx-64 conda env create --file "$env_file" --solver=libmamba
            
            # Check if the command was successful
            if [ $? -eq 0 ]; then
                # OPTIMIZATION 4: Lock the architecture so future `conda install` commands in this env don't break it
                conda run -n "$env_name" conda config --env --set subdir osx-64
                echo "  [SUCCESS] $env_name created successfully."
            else
                echo "  [ERROR] Failed to create $env_name."
            fi
        fi
    else
        echo "  [WARNING] File $env_file not found. Skipping."
    fi
    echo "------------------------------------------"
done

echo "All tasks finished."
