# Phylogenetic analysis of BA.2.86

## Installation

1. Install Nexstrain CLI following the instructions at <https://docs.nextstrain.org/projects/cli/en/stable/installation/>
2. Setup a runtime, e.g. "docker" using `nextstrain setup docker`
3. Check your installation with `nextstrain check-setup --set-default`

The easiest way to make you have all the right software versions, you can use the docker image `corneliusroemer/ba286` with `nextstrain build`.

To use that docker image, you need to add the `--docker --image=corneliusroemer/ba286` flags to the `nextstrain build` command.

Note: If you don't use the docker image, you need to make sure you have treetime v0.11.1 (or higher) installed which is at the time of writing incompatible with augur. That's why it's easiest for now to use docker and not worry.

## Getting the data

For the paper, we sampled the BA.2 background set based on an [ncov-ingest](https://www.github.com/nexstrain/ncov-ingest) curated GISAID data dump and downloaded BA.2.86 sequences from GISAID's web interface.

### BA.2 background set

Unfortunately, we cannot share the `ncov-ingest` curated data due to GISAID's data sharing policy. However, to at least partially reproduce the subsampling part of the analysis, you can use `ncov-ingest`'s "open" data, which is based on Genbank data and is freely available but less globally representative and balanced. To do so, you need to add `--config data_source=open` to the workflow invocation (or set that option in `config/config_dict.yaml`).

The other option is to use the exact BA.2 sequences from the paper. To do so, download the "Input for Augur pipeline" tarball for `EPI_SET_230907cf` from GISAID and place the archive at `data/background.tar`. Then, add `--config data_source=frozen` to the workflow invocation (or set that option in `config/config_dict.yaml`).

We generated the build using `--config data_source=gisaid`.

### BA.2.86 sequences

Besides the BA.2 data, you also need to download BA.2.86 sequences from GISAID. To reproduce the analysis in the paper with the exact BA.2.86 sequences, download `EPI_SET_230909co` and place the archive at `data/BA286.tar`. Alternatively, you can use the sequences from the paper by adding `--config data_source=frozen` to the workflow invocation (or set that option in `config/config_dict.yaml`).

Alternatively, you can download any set of BA.2.86 sequences (for example more current) and also use them for the analysis. Just download the "Input for Augur pipeline" tarball from GISAID and place the archive at `data/BA286.tar`.

## Running the workflow

Run the workflow with:

```bash
nextstrain build --docker --image corneliusroemer/ba286 .
```

You can then view the resulting tree with Auspice using:

```bash
nextstrain view auspice/BA.2.86.json
```
