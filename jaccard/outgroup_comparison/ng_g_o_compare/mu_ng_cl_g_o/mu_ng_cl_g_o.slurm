#!/usr/bin/bash
#SBATCH --partition=Tucana
#SBATCH --job-name=mu_ng_cl_g_o
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=995gb
#SBATCH --time=20-00:00:00
#SBATCH --output=%x_%j.out


# --- Go to working directory ---
cd /users/omarne/small_dataset/jaccard/outgroup_comparison/ng_g_o_compare/mu_ng_cl_g_o || exit
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
bash mu_ng_cl_g_o.sh
