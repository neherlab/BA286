regions = ["Europe", "America", "Oceania", "Asia", "Africa"]


rule build:
    input:
        "auspice/BA.2.86.json",


# rule deploy:
#     input:
#         "deploy/BA.2.86/latest",


wildcard_constraints:
    region="|".join(regions),


rule download_sequences:
    output:
        sequences="data/sequences.fasta.zst",
    shell:
        "aws s3 cp s3://nextstrain-ncov-private/aligned.fasta.zst {output.sequences}"


rule download_metadata:
    output:
        metadata="data/metadata.tsv.zst",
    shell:
        "aws s3 cp s3://nextstrain-ncov-private/metadata.tsv.zst {output.metadata}"


rule download_nextclade_dataset:
    output:
        directory("builds/nextclade_dataset"),
    shell:
        "nextclade dataset get --name='sars-cov-2' --output-dir={output}"


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


rule filter_BA286:
    input:
        metadata="builds/metadata_21L.tsv.zst",
    output:
        metadata="builds/metadata_BA286_raw.tsv",
    shell:
        """
        zstdcat {input} \
        | tsv-filter -H --str-in-fld Nextclade_pango:"BA.2.86" \
        >{output.metadata}
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
            if len(datestring) == 7:
                return datestring + "-XX"
            else:
                return datestring


        df = pd.read_csv(input.metadata, sep="\t", low_memory=False, dtype=str)
        df.loc[df.country == "Denmark", "date"] = df.loc[
            df.country == "Denmark", "date"
        ].apply(danish_daterange)
        df.date = df.date.apply(add_xx_to_incomplete_date)
        df.to_csv(output.metadata, sep="\t", index=False)


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


samples = {
    "Europe": 40,
    "America": 50,
    "Oceania": 10,
    "Asia": 80,
    "Africa": 100,
}


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


rule collect_metadata:
    input:
        regions=expand("builds/metadata_21L_{region}_subsampled.tsv", region=regions),
        BA286="builds/metadata_BA286.tsv",
        consensus_metadata="config/consensus_metadata.tsv",
    output:
        metadata="builds/metadata_sample.tsv",
        strains="builds/strains_sample.txt",
    shell:
        """
        tsv-append -H {input} >{output.metadata}
        tsv-select -H -f strain {output.metadata} >{output.strains}
        """


rule collect_sequences:
    input:
        sequences="builds/sequences_21L.fasta.zst",
        strain_names="builds/strains_sample.txt",
        consensus_sequences="config/consensus_sequences.fasta",
    output:
        sequences="builds/sequences_sample.fasta",
    shell:
        """
        mkdir -p builds/tmp
        zstdcat {input.sequences} \
        | seqkit grep -f {input.strain_names} >builds/tmp/seq.fasta
        seqkit seq -w 0 builds/tmp/seq.fasta {input.consensus_sequences} >{output.sequences}
        rm builds/tmp/seq.fasta
        """


genes = [
    "E",
    "M",
    "N",
    "ORF1a",
    "ORF1b",
    "ORF3a",
    "ORF6",
    "ORF7a",
    "ORF7b",
    "ORF8",
    "ORF9b",
    "S",
]


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
        fasta="builds/sequences_sample.fasta",
        dataset=rules.download_nextclade_dataset.output,
    output:
        alignment="builds/premask.fasta",
    shell:
        """
        nextalign run \
            --input-ref "builds/nextclade_dataset/reference.fasta" \
            --input-gene-map "builds/nextclade_dataset/genemap.gff" \
            {input.fasta} \
            --output-fasta {output.alignment} \
             2>&1
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
        dataset=rules.download_nextclade_dataset.output,
    output:
        alignment="builds/aligned.fasta",
        translations=expand("builds/translations/gene.{gene}.fasta", gene=genes),
    params:
        template_string=lambda w: "builds/translations/gene.{gene}.fasta",
    shell:
        """
        nextalign run \
            --input-ref "builds/nextclade_dataset/reference.fasta" \
            --input-gene-map "builds/nextclade_dataset/genemap.gff" \
            --genes E,M,N,ORF1a,ORF1b,ORF3a,ORF6,ORF7a,ORF7b,ORF8,ORF9b,S \
            {input.fasta} \
            --output-translations {params.template_string} \
            --output-fasta {output.alignment} \
             2>&1
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
            --root "BA.2" \
            --output {output.tree}
        """


rule refine:
    input:
        tree="builds/tree_fixed.nwk",
        alignment="builds/aligned.fasta",
        metadata="builds/metadata.tsv",
    output:
        tree="builds/tree.nwk",
        node_data="builds/branch_lengths.json",
    shell:
        """
        augur refine \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --divergence-units mutations \
            --output-tree {output.tree} \
            --metadata {input.metadata} \
            --output-node-data {output.node_data} \
            --timetree \
            --coalescent skyline \
            --keep-polytomies \
            --date-inference marginal \
            --root "BA.2" \
            --date-confidence \
            --no-covariance
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


def _get_node_data_by_wildcards(wildcards):
    """Return a list of node data files to include for a given build's wildcards."""
    # Define inputs shared by all builds.
    wildcards_dict = dict(wildcards)
    inputs = [
        rules.refine.output.node_data,
    ]
    inputs = [input_file.format(**wildcards_dict) for input_file in inputs]
    return inputs


rule export:
    input:
        tree="builds/tree.nwk",
        node_data=_get_node_data_by_wildcards,
        ancestral="builds/muts.json",
        auspice_config="config/auspice_config.json",
        lat_longs=rules.download_lat_longs.output,
        metadata="builds/metadata.tsv",
        description="config/description.md",
    output:
        auspice_json="auspice/BA.2.86.json",
        root_json="auspice/BA.2.86_root-sequence.json",
    shell:
        """
        augur export v2 \
            --tree {input.tree} \
            --node-data {input.node_data} {input.ancestral} \
            --include-root-sequence \
            --auspice-config {input.auspice_config} \
            --lat-longs {input.lat_longs} \
            --output {output.auspice_json} \
            --metadata {input.metadata} \
            --description {input.description}
        """


# rule deploy_single:
#     input:
#         "auspice/ncov_{build}.json",
#         "auspice/ncov_{build}_root-sequence.json",
#     output:
#         "deploy/{build}/latest",
#     shell:
#         """
#         nextstrain remote upload nextstrain.org/groups/neherlab/ncov/BA.2.86 {input} 2>&1 && touch {output}
#         """
