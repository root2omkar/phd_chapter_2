#!/usr/bin/bash
#SBATCH --partition=Tucana
#SBATCH --job-name=ma_mu_gi
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=995gb
#SBATCH --time=10-00:00:00
#SBATCH --output=%x_%j.out


# --- Go to working directory ---
cd /users/omarne/small_dataset/jaccard/ingroup_comparison/g_g_i_compare/ma_g_mu_g_i || exit
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
bash ma_g_mu_g_i.sh
