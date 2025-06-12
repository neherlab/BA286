data_source = config["data_source"]
regions = config["regions"]
samples = config["samples"]
genes = config["genes"]

if data_source == "frozen":
    use_frozen_background = True
else:
    use_frozen_background = False

if data_source == "open":
    use_open = True
    prefix = "https://data.nextstrain.org/files/ncov/open/"
else:
    use_open = False
    prefix = "s3://nextstrain-ncov-private/"


rule build:
    input:
        "auspice/BA.3.2.json",


wildcard_constraints:
    region="|".join(regions),


rule download_sequences:
    output:
        sequences="data/sequences.fasta.zst",
    params:
        url=prefix + "aligned.fasta.zst",
        open=use_open,
    shell:
        """
        if [ {params.open} == True ]; then
            curl {params.url} -o {output.sequences}
        else
            aws s3 cp {params.url} {output.sequences}
        fi
        """


rule download_metadata:
    output:
        metadata="data/metadata.tsv.zst",
    params:
        url=prefix + "metadata.tsv.zst",
        open=use_open,
    shell:
        """
        if [ {params.open} == True ]; then
            curl {params.url} -o {output.metadata}
        else
            aws s3 cp {params.url} {output.metadata}
        fi
        """


rule download_lat_longs:
    output:
        "builds/lat_longs.tsv",
    params:
        url="https://raw.githubusercontent.com/nextstrain/ncov/master/defaults/lat_longs.tsv",
    shell:
        """
        curl {params.url} | \
        sed "s/North Rhine Westphalia/North Rhine-Westphalia/g" | \
        sed "s/Baden-Wuerttemberg/Baden-Wurttemberg/g" \
        > {output}
        """


rule filter_21L:
    input:
        metadata="data/metadata.tsv.zst",
    output:
        metadata="builds/metadata_21L.tsv.zst",
        strains="builds/strains_21L.txt",
    shell:
        """
        zstdcat {input} \
        | tsv-filter -H --str-eq clade_nextstrain:21L \
        | zstd >{output.metadata}
        zstdcat {output.metadata} \
        | tsv-select -H -f strain \
        >{output.strains}
        """


rule filter_21L_sequences:
    input:
        sequences="data/sequences.fasta.zst",
        strains="builds/strains_21L.txt",
    output:
        sequences="builds/sequences_21L.fasta.zst",
    shell:
        """
        zstdcat {input} \
        | seqkit grep -f {input.strains} \
        | zstd >{output.sequences}
        """


rule filter_21L_region:
    input:
        metadata="builds/metadata_21L.tsv.zst",
    output:
        metadata="builds/metadata_21L_{region}.tsv.zst",
    shell:
        """
        zstdcat {input} \
        | tsv-filter -H --str-in-fld "region:{wildcards.region}" \
        | zstd >{output.metadata}
        """


rule subsample_21L_region:
    input:
        metadata="builds/metadata_21L_{region}.tsv.zst",
    output:
        metadata="builds/metadata_21L_{region}_subsampled.tsv",
    params:
        number=lambda w: samples[w.region],
    shell:
        """
        augur filter \
            --metadata {input.metadata} \
            --min-date 2021-11-01 \
            --max-date 2022-06-01 \
            --exclude-ambiguous-dates-by any \
            --exclude-where "reversion_mutations!=0" "QC_stop_codons!=good" "QC_frame_shifts!=good" "QC_snp_clusters!=good" "QC_rare_mutations!=good" "QC_missing_data!=good" \
            --subsample-max-sequences {params.number} \
            --subsample-seed 42 \
            --output-metadata {output.metadata}
        """


rule get_tarball:
    """
    Copy tarball from ~/Downloads to data directory
    Get the most recent tarball (by date)
    """
    output:
        tarball="data/gisaid.tar",
    params:
        lastest_tarball_cmd='find /Users/corneliusromer/Downloads/ -maxdepth 1 -name "gisaid_auspice_input_hcov-19_202*.tar" | sed "s/.\{4\}$//" | sort | tail -1 | sed "s/$/.tar/"',
    shell:
        """
        latest_tarball=$(find /Users/corneliusromer/Downloads/ -maxdepth 1 -name "gisaid_auspice_input_hcov-19_202*.tar" | sed "s/....$//" | sort | tail -1 | sed "s/$/.tar/")
        echo "$latest_tarball"
        cp "$latest_tarball" {output.tarball}
        """


rule unpack_BA286:
    """
    Need to download the tarball from GISAID and place it in the data directory.
    """
    input:
        "data/gisaid.tar",
    output:
        metadata="builds/metadata_BA286_raw.tsv",
        sequences="builds/sequences_BA286.fasta",
    params:
        cache="data/cache",
    shell:
        """
        mkdir -p builds/tmp
        tar -xvf {input} -C builds/tmp
        cp {params.cache}/* builds/tmp
        tsv-uniq -H -f strain builds/tmp/*.tsv >{output.metadata}
        cat builds/tmp/*.fasta | seqkit rmdup -n >{output.sequences}
        # rm -rf builds/tmp
        cp {output.metadata} {params.cache}
        cp {output.sequences} {params.cache}
        echo $(seqkit seq -n {output.sequences} | wc -l) BA.3.2 sequences
        echo $(tail -n +2 {output.metadata} | wc -l) BA.3.2 metadata entries
        """


rule unpack_background:
    """
    Frozen background for reproducibility
    """
    input:
        tarball="data/background.tar",
    output:
        metadata="builds/metadata_background.tsv",
        sequences="builds/sequences_background.fasta",
    shell:
        """
        mkdir -p builds/tmpbackground
        tar -xvf {input.tarball} -C builds/tmpbackground
        mv builds/tmpbackground/*.metadata.tsv {output.metadata}
        mv builds/tmpbackground/*.fasta {output.sequences}
        rm -rf builds/tmpbackground
        """


rule postprocess_dates:
    """
    - Add XX to incomplete dates
    - Convert Danish dates to date ranges
    """
    input:
        metadata="builds/metadata_BA286_raw.tsv",
    output:
        metadata="builds/metadata_BA286.tsv",
    run:
        import pandas as pd
        from treetime.utils import numeric_date
        from datetime import datetime as dt
        from datetime import timedelta


        def danish_daterange(datestring):
            start_date = dt.strptime(datestring, "%Y-%m-%d")
            end_date = start_date + timedelta(days=7)
            return f"[{numeric_date(start_date)}:{numeric_date(end_date)}]"


        def add_xx_to_incomplete_date(datestring):
            if len(datestring) == 4:
                return datestring + "-XX-XX"
            elif len(datestring) == 7:
                return datestring + "-XX"
            else:
                return datestring


        df = pd.read_csv(input.metadata, sep="\t", low_memory=False, dtype=str)
        df.loc[df.country == "Denmark", "date"] = df.loc[
            df.country == "Denmark", "date"
        ].apply(danish_daterange)
        df.date = df.date.apply(add_xx_to_incomplete_date)
        df.to_csv(output.metadata, sep="\t", index=False)


rule collect_metadata:
    input:
        background=lambda w: "builds/metadata_background.tsv"
        if use_frozen_background
        else expand("builds/metadata_21L_{region}_subsampled.tsv", region=regions),
        BA286="builds/metadata_BA286.tsv",
        consensus_metadata="config/consensus_metadata.tsv",
    output:
        metadata="builds/metadata_sample_raw.tsv",
        strains="builds/strains_sample.txt",
    params:
        fields="strain,gisaid_epi_isl,date,region,country,division,location,host,age,sex,originating_lab,submitting_lab,authors,url,title,paper_url,date_submitted",
    shell:
        """
        rm -rf builds/tmpcollect
        mkdir -p builds/tmpcollect
        for file in {input}; do
            tsv-select -H -f {params.fields} $file >builds/tmpcollect/$(basename $file)
        done
        tsv-append -H builds/tmpcollect/*.tsv >{output.metadata}
        tsv-select -H -f strain {output.metadata} >{output.strains}
        """


rule rename_southern_africa:
    """
    Add a new column "renamed_region" to the metadata
    This column is like "region" but if country is South Africa, region is "Southern Africa"
    """
    input:
        metadata="builds/metadata_sample_raw.tsv",
    output:
        metadata="builds/metadata_sample.tsv",
    run:
        import pandas as pd

        df = pd.read_csv(input.metadata, sep="\t", low_memory=False, dtype=str)
        df["renamed_region"] = df.region
        df.loc[df.renamed_region == "Africa", "renamed_region"] = "Other Africa"
        df.loc[
            df["country"].isin(["South Africa", "Eswatini", "Namibia"]),
            "renamed_region",
            ] = "Southern Africa"
        df.to_csv(output.metadata, sep="\t", index=False)


rule get_sampled_sequences:
    input:
        sequences="builds/sequences_21L.fasta.zst",
        strain_names="builds/strains_sample.txt",
    output:
        sequences="builds/sequences_sampled.fasta",
    shell:
        """
        zstdcat {input.sequences} \
        | seqkit grep -f {input.strain_names} >{output.sequences}
        """


rule collect_sequences:
    input:
        background=lambda w: "builds/sequences_background.fasta"
        if use_frozen_background
        else "builds/sequences_sampled.fasta",
        consensus_sequences="config/consensus_sequences.fasta",
        BA286="builds/sequences_BA286.fasta",
    output:
        sequences="builds/sequences_sample.fasta",
    shell:
        """
        seqkit seq -w 0 {input} >{output.sequences}
        """


rule epi_sets:
    input:
        ba286="builds/metadata_BA286.tsv",
        total="builds/metadata_sample_raw.tsv",
    output:
        ba286="builds/epi_sets/BA286.txt",
        total="builds/epi_sets/total.txt",
    shell:
        """
        tsv-select -H -f gisaid_epi_isl {input.ba286} | tail -n +2 >{output.ba286}
        tsv-select -H -f gisaid_epi_isl {input.total} | tail -n +2 >{output.total}
        """


rule exclude_outliers:
    input:
        sequences="builds/sequences_sample.fasta",
        metadata="builds/metadata_sample.tsv",
        exclude="config/exclude.txt",
    output:
        sequences="builds/sequences.fasta",
        metadata="builds/metadata.tsv",
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --exclude {input.exclude} \
            --output {output.sequences} \
            --output-metadata {output.metadata}
        """


rule nextalign_before_mask:
    input:
        fasta="builds/sequences.fasta",
        reference="config/reference.fasta",
        genemap="config/genemap.gff",
    output:
        alignment="builds/premask.fasta",
    shell:
        """
        nextalign run \
            --input-ref {input.reference} \
            --input-gene-map {input.genemap} \
            --output-fasta {output.alignment} \
            -- {input.fasta}
        """


rule mask_specific:
    input:
        alignment="builds/premask.fasta",
        mask_file="config/mask.tsv",  # I assume you want to input the mask file in this format.
    output:
        alignment="builds/masked_specific.fasta",
    shell:
        """
        python3 scripts/mask_specific.py \
            --input-alignment {input.alignment} \
            --mask-file {input.mask_file} \
            --output-alignment {output.alignment}
        """


rule mask:
    input:
        alignment="builds/masked_specific.fasta",
        mask_sites="config/mask_sites.txt",
    output:
        alignment="builds/masked.fasta",
    shell:
        """
        python3 scripts/mask-alignment.py \
            --alignment {input.alignment} \
            --mask-site-file {input.mask_sites} \
            --mask-from-beginning 100 \
            --mask-from-end 200 \
            --mask-terminal-gaps \
            --mask-all-gaps \
            --output {output.alignment}
        """


rule nextalign_after_mask:
    input:
        fasta="builds/masked.fasta",
        reference="config/reference.fasta",
        genemap="config/genemap.gff",
    output:
        alignment="builds/aligned.fasta",
        translations=expand("builds/translations/gene.{gene}.fasta", gene=genes),
    params:
        template_string=lambda w: "builds/translations/gene.{gene}.fasta",
    shell:
        """
        nextalign run \
            --input-ref {input.reference} \
            --input-gene-map {input.genemap} \
            --genes E,M,N,ORF1a,ORF1b,ORF3a,ORF6,ORF7a,ORF7b,ORF8,ORF9b,S \
            --output-translations {params.template_string} \
            --output-fasta {output.alignment} \
            -- {input.fasta}
        """


rule tree:
    input:
        alignment="builds/aligned.fasta",
        exclude_sites="config/tree_exclude_sites.txt",
    output:
        tree="builds/tree_raw.nwk",
    params:
        args="'-ninit 10 -n 4 -czb -T AUTO'",
    shell:
        """
        augur tree \
            --exclude-sites {input.exclude_sites} \
            --alignment {input.alignment} \
            --tree-builder-args {params.args} \
            --output {output.tree}
        """


rule fix_iqtree:
    input:
        tree="builds/tree_raw.nwk",
        alignment="builds/aligned.fasta",
    output:
        tree="builds/tree_fixed.nwk",
    shell:
        """
        python3 scripts/fix_tree.py \
            --alignment {input.alignment} \
            --input-tree {input.tree} \
            --root "BA.3" \
            --output {output.tree}
        """


rule create_figure:
    input:
        metadata="builds/metadata.tsv",
        tree="builds/tree_fixed.nwk",
        alignment="builds/aligned.fasta",
    output:
        figure="figures/timetree.pdf",
        node_data="builds/branch_lengths.json",
        tree="builds/tree.nwk",
    shell:
        """
        python3 scripts/timetree_inference.py \
            --metadata {input.metadata} \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --figure {output.figure} \
            --output-node-data {output.node_data} \
            --output-tree {output.tree}
        """


rule ancestral:
    input:
        tree="builds/tree.nwk",
        alignment="builds/aligned.fasta",
        annotation="config/annotation.gb",
        translations=expand("builds/translations/gene.{gene}.fasta", gene=genes),
    output:
        node_data="builds/muts.json",
    params:
        inference="joint",
        translations="builds/translations/gene.%GENE.fasta",
    shell:
        """
        augur ancestral \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --output-node-data {output.node_data} \
            --inference {params.inference} \
            --genes E M N ORF1a ORF1b ORF3a ORF6 ORF7a ORF7b ORF8 ORF9b S \
            --annotation {input.annotation} \
            --translations {params.translations} \
            2>&1 | tee {log}
        """


rule colors:
    input:
        ordering="config/color_ordering.tsv",
        color_schemes="config/color_schemes.tsv",
        metadata="builds/metadata.tsv",
    output:
        colors="builds/colors.tsv",
    shell:
        """
        python3 scripts/assign-colors.py \
            --ordering {input.ordering} \
            --color-schemes {input.color_schemes} \
            --output {output.colors} \
            --metadata {input.metadata} 2>&1
        """


rule export:
    input:
        tree="builds/tree.nwk",
        node_data="builds/branch_lengths.json",
        ancestral="builds/muts.json",
        auspice_config="config/auspice_config.json",
        lat_longs=rules.download_lat_longs.output,
        metadata="builds/metadata.tsv",
        colors="builds/colors.tsv",
        description="config/description.md",
    output:
        auspice_json="builds/BA.3.2.json",
        root_json="builds/BA.3.2_root-sequence.json",
    shell:
        """
        augur export v2 \
            --tree {input.tree} \
            --node-data {input.node_data} {input.ancestral} \
            --colors {input.colors} \
            --include-root-sequence \
            --auspice-config {input.auspice_config} \
            --lat-longs {input.lat_longs} \
            --output {output.auspice_json} \
            --metadata {input.metadata} \
            --description {input.description}
        """


rule postprocessing:
    input:
        auspice="builds/BA.3.2.json",
        root="builds/BA.3.2_root-sequence.json",
    output:
        auspice="auspice/BA.3.2.json",
        root="auspice/BA.3.2_root-sequence.json",
    shell:
        """
        python3 scripts/postprocess.py
        cp {input.root} {output.root}
        """


rule deploy:
    input:
        "auspice/BA.3.2.json",
        "auspice/BA.3.2_root-sequence.json",
    shell:
        """
        nextstrain remote upload nextstrain.org/groups/neherlab/ncov/BA.3.2 {input} 2>&1
        """
