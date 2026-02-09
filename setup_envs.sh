#!/bin/bash

# Define the list of environment files
# Add or remove files here as needed
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
echo "Starting Conda Environment Setup"
echo "=========================================="

# check if conda is installed
if ! command -v conda &> /dev/null; then
    echo "Error: Conda is not installed or not in your PATH."
    exit 1
fi

# Loop through each file
for env_file in "${ENV_FILES[@]}"; do
    if [ -f "$env_file" ]; then
        # Extract environment name from the YAML file to check if it exists
        env_name=$(grep "name:" "$env_file" | head -n 1 | cut -d ' ' -f 2)
        
        echo "Processing $env_file (Env Name: $env_name)..."

        # Check if environment already exists
        if conda info --envs | grep -q "^$env_name "; then
            echo "  [SKIP] Environment '$env_name' already exists."
        else
            echo "  [CREATE] Creating environment '$env_name'..."
            conda env create --file "$env_file" --solver=libmamba -y
            
            # Check if the command was successful
            if [ $? -eq 0 ]; then
                echo "  [SUCCESS] $env_name created successfully."
            else
                echo "  [ERROR] Failed to create $env_name."
                # Optional: exit 1 # Uncomment to stop script on first failure
            fi
        fi
    else
        echo "  [WARNING] File $env_file not found. Skipping."
    fi
    echo "------------------------------------------"
done

echo "All tasks finished.":
