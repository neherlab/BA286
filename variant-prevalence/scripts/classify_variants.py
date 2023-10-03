#%%
"""
Reads in TSV with each row being singl observation
- Classifies each sample as belonging to a variant
- Aggregates into weekly counts
"""
import numpy as np
import pandas as pd
from pango_aliasor.aliasor import Aliasor
aliasor = Aliasor()
#%%
df = pd.read_csv('results/metadata-south-africa.tsv', sep='\t')
df

#%%
#
def safe_float(x):
    try:
        return float(x)
    except:
        return np.nan
# Remove samples that have NAN in Nextclade_pango or clade_nextstrain
df.clock_deviation = df["clock_deviation"].apply(safe_float)
df.dropna(subset=['Nextclade_pango', 'clade_nextstrain', "clock_deviation"], inplace=True)
df

#%%
# Use clock deviation: needs to be >= -10 and <= 25
df = df[(df.clock_deviation >= -10) & (df.clock_deviation <= 20)]

#%%
# Overwrite NextcladePango with unaliased version
df['Nextclade_pango'] = df['Nextclade_pango'].apply(lambda x: aliasor.uncompress(x))
df
#%%
# Variant mapping
VARIANT_MAP = {
    "Beta": df.clade_nextstrain.isin(["20H"]),
    "Alpha": df.clade_nextstrain.isin(["20I"]),
    "Delta": df.clade_nextstrain.isin(["21A", "21I", "21J"]),
    "BA.1": df.clade_nextstrain.isin(["21K", "21M"]),
    "BA.2": df.clade_nextstrain.isin(["21L", "22C"]),
    "BA.2.75": df.clade_nextstrain.isin(["22D", "23C"]),
    "BA.4/5": df.clade_nextstrain.isin(["22A", "22B", "22E"]),
    "XBB": df.Nextclade_pango.str.startswith("XBB"),
    "BA.2.86": df.Nextclade_pango.str.startswith("B.1.1.529.2.86"),
}
df["variant"] = "Other"
for variant, mask in VARIANT_MAP.items():
    df.loc[mask, "variant"] = variant
df
# %%

# Round date down to Monday of week, output format as YYYY-MM-DD
df["week"] = pd.to_datetime(df.date).dt.to_period("W-MON").dt.strftime("%Y-%m-%d")
df
# %%
# Aggregate counts by week in wide format
wide = df.groupby(["week", "variant"]).size().unstack().fillna(0).astype(int)
#%%
# Add total column
wide["total"] = wide.sum(axis=1)

# Add name of most common variant column
wide["most_common"] = wide.drop(columns=["total"]).idxmax(axis=1)

# Add count of most common variant column
wide["most_common_count"] = wide.drop(columns=["total"]).max(axis=1)

# Add percentag of most common variant column
wide["most_common_pct"] = 100*wide["most_common_count"] / wide["total"]

# Add percentage columns
from itertools import chain
for variant in chain(VARIANT_MAP.keys(),["Other"]):
    wide[f"{variant}_pct"] = wide[variant] / wide["total"] * 100


# %%
COLUMNS=[
    "total",
    "most_common",
    "most_common_count",
    "most_common_pct",
    "Other",
    "Beta",
    "Alpha",
    "Delta",
    "BA.1",
    "BA.2",
    "BA.4/5",
    "BA.2.75",
    "XBB",
    "BA.2.86",
    "Other_pct",
    "Beta_pct",
    "Alpha_pct",
    "Delta_pct",
    "BA.1_pct",
    "BA.2_pct",
    "BA.4/5_pct",
    "BA.2.75_pct",
    "XBB_pct",
    "BA.2.86_pct",
]
# wide.to_csv("results/weekly-variant-counts-south-africa.tsv", sep="\t", columns=COLUMNS)
# format floats with 2 decimal places
wide.to_csv("results/weekly-variant-counts-south-africa.tsv", sep="\t", columns=COLUMNS, float_format="%.1f")

# %%
# Make plot of all the pct columns
# wide[[col for col in wide.columns if (col.endswith("_pct") and col != "most_common_pct")]].plot.line(figsize=(12, 6), title="SARS-CoV-2 variant proportions in South Africa")

# %%
