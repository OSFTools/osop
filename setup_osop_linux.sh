# ============================================================
# Script Title: osop download and run
#
# 1. Dr. Amal Saber ElSehwagy
#    Email: amal.saber.1192@gmail.com
#
# 2. Dr. Marwa S. Mohamed
#    Email: marwa.samy.ema@outlook.com
#
# Year: 2026
# ============================================================

#!/bin/bash

# Clone OSOP
git clone https://github.com/OSFTools/osop.git

# Download and install Miniconda
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b -p "$HOME/miniconda3"

# Initialize Conda
"$HOME/miniconda3/bin/conda" init bash
source "$HOME/.bashrc"

# Accept Anaconda Terms of Service
"$HOME/miniconda3/bin/conda" tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main
"$HOME/miniconda3/bin/conda" tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r

# Create OSOP environment
cd "$HOME/osop/envs"
"$HOME/miniconda3/bin/conda" create --file osop-linux-spec-file.yml --name osop -y

# Activate environment
source "$HOME/.bashrc"
conda activate osop
