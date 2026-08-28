#!/usr/bin/env python3

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats

# -----------------------------
# LOAD DATA
# -----------------------------
pos_df = pd.read_csv("jac_i_align_pos.csv", index_col=0)
states_df = pd.read_csv("jac_i_align_states.csv", index_col=0)

# Convert to numeric
pos_df = pos_df.apply(pd.to_numeric, errors='coerce')
states_df = states_df.apply(pd.to_numeric, errors='coerce')


# -----------------------------
# HEATMAP FUNCTION (LOWER TRIANGLE)
# -----------------------------
def plot_lower_triangle_heatmap(df, title, filename):
    mask = np.triu(np.ones_like(df, dtype=bool))

    plt.figure(figsize=(20, 15))

    # Gray colormap with light gray instead of pure white for 0
    cmap = sns.light_palette("gray", as_cmap=True)

    sns.heatmap(
        df,
        mask=mask,
        cmap=cmap,
        annot=True,
        fmt=".2f",
        annot_kws={"size": 22},
        square=True,
        linewidths=0.5,
        linecolor="black",
        cbar_kws={"label": "Jaccard similarity"}
    )

    plt.title(title, fontsize=22)
    plt.xticks(fontsize=22)
    plt.yticks(fontsize=22)

    # Increase colorbar font size
    cbar = plt.gcf().axes[-1]
    cbar.set_ylabel("Jaccard similarity")
    cbar.tick_params(labelsize=22)

    # Add border around the heatmap
    ax = plt.gca()
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("black")
        spine.set_linewidth(1.5)

    plt.tight_layout()
    plt.savefig(filename, format="pdf")
    plt.close()


# Generate heatmaps
plot_lower_triangle_heatmap(pos_df,
                           "Positions", 
                           "jaccard_heatmap_ingroup_positions.pdf")

plot_lower_triangle_heatmap(states_df,
                           "States", 
                           "jaccard_heatmap_ingroup_states.pdf")


# -----------------------------
# EXTRACT CLEAN VALUES
# -----------------------------
def get_values(df):
    mask = np.tril(np.ones_like(df, dtype=bool), k=-1)
    values = df.where(mask).stack().values

    # Remove NaNs
    values = values[~np.isnan(values)]
    return values


pos_values = get_values(pos_df)
states_values = get_values(states_df)


# -----------------------------
# FORMAT FUNCTION
# -----------------------------
def fmt(x):
    if isinstance(x, (float, np.floating)):
        if np.isnan(x):
            return "NaN"
        return f"{x:.6f}"
    return str(x)


# -----------------------------
# LOGGING SETUP
# -----------------------------
results = []

def log(msg):
    print(msg)
    results.append(msg)


log("===== STATISTICAL ANALYSIS REPORT =====\n")
log(f"Positions N = {len(pos_values)}")
log(f"States N = {len(states_values)}")


# -----------------------------
# NORMALITY TESTS
# -----------------------------
log("\n--- NORMALITY TESTS ---")

if len(pos_values) >= 3:
    shapiro_pos = stats.shapiro(pos_values)
    log(f"Shapiro (Positions): statistic={fmt(shapiro_pos.statistic)}, p={fmt(shapiro_pos.pvalue)}")
else:
    log("Shapiro (Positions): Not enough data")

if len(states_values) >= 3:
    shapiro_states = stats.shapiro(states_values)
    log(f"Shapiro (States): statistic={fmt(shapiro_states.statistic)}, p={fmt(shapiro_states.pvalue)}")
else:
    log("Shapiro (States): Not enough data")

# -----------------------------
# Kolgomorov-Smirnov test
# -----------------------------
ks_pos = stats.kstest(pos_values, 'norm', args=(np.mean(pos_values), np.std(pos_values)))
ks_states = stats.kstest(states_values, 'norm', args=(np.mean(states_values), np.std(states_values)))

log(f"KS (Positions): statistic={fmt(ks_pos.statistic)}, p={fmt(ks_pos.pvalue)}")
log(f"KS (States): statistic={fmt(ks_states.statistic)}, p={fmt(ks_states.pvalue)}")


# -----------------------------
# VARIANCE TESTS
# -----------------------------
log("\n--- VARIANCE TESTS ---")

levene_test = stats.levene(pos_values, states_values)
bartlett_test = stats.bartlett(pos_values, states_values)

log(f"Levene: statistic={fmt(levene_test.statistic)}, p={fmt(levene_test.pvalue)}")
log(f"Bartlett: statistic={fmt(bartlett_test.statistic)}, p={fmt(bartlett_test.pvalue)}")


# -----------------------------
# PARAMETRIC TEST
# -----------------------------
log("\n--- PARAMETRIC TEST ---")

t_test = stats.ttest_ind(pos_values, states_values, equal_var=False)
log(f"Welch T-test: statistic={fmt(t_test.statistic)}, p={fmt(t_test.pvalue)}")


# -----------------------------
# NON-PARAMETRIC TESTS
# -----------------------------
log("\n--- NON-PARAMETRIC TESTS ---")

mann_whitney = stats.mannwhitneyu(pos_values, states_values, alternative='two-sided')
log(f"Mann–Whitney U: statistic={fmt(mann_whitney.statistic)}, p={fmt(mann_whitney.pvalue)}")

if len(pos_values) == len(states_values) and len(pos_values) > 0:
    wilcoxon = stats.wilcoxon(pos_values, states_values)
    log(f"Wilcoxon: statistic={fmt(wilcoxon.statistic)}, p={fmt(wilcoxon.pvalue)}")
else:
    log("Wilcoxon: Skipped (unequal lengths or insufficient data)")


# -----------------------------
# DISTRIBUTION TESTS
# -----------------------------
log("\n--- DISTRIBUTION TESTS ---")

ks_2sample = stats.ks_2samp(pos_values, states_values)
log(f"KS 2-sample: statistic={fmt(ks_2sample.statistic)}, p={fmt(ks_2sample.pvalue)}")

anderson = stats.anderson_ksamp([pos_values, states_values])
log(f"Anderson–Darling: statistic={fmt(anderson.statistic)}, p={fmt(anderson.pvalue)}")


# -----------------------------
# EFFECT SIZE
# -----------------------------
log("\n--- EFFECT SIZE ---")

def cohens_d(x, y):
    if len(x) < 2 or len(y) < 2:
        return np.nan
    nx, ny = len(x), len(y)
    pooled_std = np.sqrt(((nx-1)*np.var(x) + (ny-1)*np.var(y)) / (nx+ny-2))
    if pooled_std == 0:
        return np.nan
    return (np.mean(x) - np.mean(y)) / pooled_std

d = cohens_d(pos_values, states_values)

u_stat = mann_whitney.statistic
n1, n2 = len(pos_values), len(states_values)

if n1 > 0 and n2 > 0:
    rbc = 1 - (2*u_stat) / (n1*n2)
else:
    rbc = np.nan

log(f"Cohen's d: {fmt(d)}")
log(f"Rank-biserial correlation: {fmt(rbc)}")


# -----------------------------
# SUMMARY STATISTICS
# -----------------------------
log("\n--- SUMMARY STATISTICS ---")

def summary(x):
    if len(x) == 0:
        return "No data"
    return {
        "mean": fmt(np.mean(x)),
        "median": fmt(np.median(x)),
        "std": fmt(np.std(x)),
        "min": fmt(np.min(x)),
        "max": fmt(np.max(x))
    }

log(f"Positions: {summary(pos_values)}")
log(f"States: {summary(states_values)}")


# -----------------------------
# SAVE RESULTS
# -----------------------------
with open("jac_ingroup_stats.txt", "w") as f:
    for line in results:
        f.write(line + "\n")

print("\n Analysis complete!")
print("📁 Files generated:")
print(" - heatmap_positions.pdf")
print(" - heatmap_states.pdf")
print(" - jac_ingroup_stats.txt")

