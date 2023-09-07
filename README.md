# Phylogenetic analysis of BA.2.86

## Installation

- Nextstrain CLI with docker backend
- TODO: Which docker image?

Requires treetime master (or at least v0.11.1 which is not yet released)

## Usage

Run the workflow with:

```bash
nextstrain build .
```

Result is viewable via:

```bash
nextstrain view auspice/BA.2.86.json
```

Set `use_frozen_background` to `false` to generate a fresh sample of background BA.2 sequences.
Otherwise, you need to put background sequences into `data/background.tar`.

## Data

Data used is from GISAID, curated using ncov-ingest

As the data is not public, we cannot share it here. However, it is possible to run the analysis using Nextstrain-curated Genbank data (TODO: Howto)

TODO: Allow reproducibility from selected GISAID sequences (using EPI set)

## EPI_sets

Background set: `EPI_SET_230907cf`, download "Input for Augur" and place in `data/background.tar`
BA.2.86 sequences: download "Input for Augur pipeline" and place into `data/gisaid.tar`