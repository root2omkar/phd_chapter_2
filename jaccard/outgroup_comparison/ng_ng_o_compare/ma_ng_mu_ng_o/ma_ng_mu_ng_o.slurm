#!/usr/bin/bash
#SBATCH --partition=Orion
#SBATCH --job-name=ma_mu_ng_o
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=3995gb
#SBATCH --time=30-00:00:00
#SBATCH --output=%x_%j.out


# --- Go to working directory ---
cd /users/omarne/small_dataset/jaccard/outgroup_comparison/ng_ng_o_compare/ma_ng_mu_ng_o || exit
pwd

# --- Initialize Conda ---
# Source the conda setup script
source /apps/pkg/anaconda3/2023.09/etc/profile.d/conda.sh

# Activate your environment
conda activate omkar

# Check Python
which python3
python3 --version

# Test imports
python -c "import parsl; import Bio; print('All OK')"

# Run jaccard script
bash ma_ng_mu_ng_o.sh


