#!/bin/bash

# Check for ARM64 architecture
ARCH=$(uname -m)
if [ "$ARCH" != "aarch64" ] && [ "$ARCH" != "arm64" ]; then
    echo "This script is intended for ARM64 systems."
    exit 1
fi

# Create a Conda environment
echo "Creating the 'puppy' Conda environment..."
conda create -y -n puppy python=3.10.6 pyarrow=14.0.1

# Install dependencies
echo "Installing dependencies in the Conda environment..."
conda run -n puppy conda install -y -c bioconda -c conda-forge mmseqs2 pandas=1.5 biopython dask matplotlib
conda run -n puppy pip install colorama seaborn primer3-py

# Install mmseqs2 for ARM64
echo "Installing mmseqs2..."
MMSEQS2_URL="https://mmseqs.com/latest/mmseqs-linux-arm64.tar.gz"
wget -qO- $MMSEQS2_URL | tar xvz -C /usr/local/bin

# Install mmseqs2 for ARM64
echo "Installing mmseqs2..."
MMSEQS2_URL="https://mmseqs.com/latest/mmseqs-linux-arm64.tar.gz"
wget -qO- $MMSEQS2_URL | tar xvz -C /usr/local/bin

# Add mmseqs2 and Puppy scripts to PATH
echo "Adding mmseqs2 and Puppy scripts to PATH..."
SCRIPTS_DIR="$(pwd)/scripts"
echo "export PATH=\$PATH:/usr/local/bin:$SCRIPTS_DIR" >> ~/.bashrc
source ~/.bashrc

echo "Installation complete! Please restart your terminal or run 'source ~/.bashrc' to apply changes."
