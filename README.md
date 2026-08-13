<div align="center">
  <h2>🧩 <code>immundapter</code> — Run immune-repertoire tools from R</h2>
</div>

<div align="center">
  <a href="https://github.com/immunomind">
    <img alt="Ecosystem: ImmunoMind"
         src="https://img.shields.io/badge/ecosystem-ImmunoMind-orange?style=flat-square">
  </a>
  <a href="https://github.com/immunomind/immundapter/actions/workflows/R-CMD-check.yaml">
    <img alt="R CMD check"
         src="https://github.com/immunomind/immundapter/actions/workflows/R-CMD-check.yaml/badge.svg">
  </a>
  <a href="https://github.com/immunomind/immundapter/issues">
    <img alt="GitHub Issues"
         src="https://img.shields.io/github/issues/immunomind/immundapter?style=flat-square">
  </a>
</div>

<p align="center">
  <a href="https://immunomind.github.io/docs/">Documentation</a>
  |
  <a href="https://github.com/immunomind/immundapter/issues">Issues</a>
  |
  <a href="https://github.com/immunomind">ImmunoMind ecosystem</a>
</p>

`immundapter` connects [`immundata`](https://github.com/immunomind/immundata)
to external tools for immune-repertoire analysis. It prepares the input, runs
the tool, and gives access to the result tables and files.

> [!WARNING]
> `immundapter` is under development. Functions may change before version 1.0.

## Installation

```r
install.packages("pak")
pak::pkg_install("immunomind/immundapter")
```

Each external tool must also be installed. Use `tool$help()` to print its
installation instructions and `tool$check()` to check that it is ready.

## The common workflow

Analysis tools use four commands:

```r
config <- conf_tcrdist3()       # Select the analysis settings.
tool <- tool_tcrdist3(config)   # Create the tool.
run <- tool$run(input)          # Run one analysis.
results <- run$results()        # Read its results.
```

`tool$run()` returns information about the completed run, including its output
directory. It does not return a biological result table directly. Use
`run$results()` to access the tables and files.

Large tables stay on disk. Use `dplyr::collect()` when you want to read one
into R:

```r
clusters <- dplyr::collect(results$clusters)
```

## Available tools

| Tool | Input passed to `tool$run()` | Main result |
| --- | --- | --- |
| tcrdist3 | `ImmunData` or a CSV/TSV/TXT receptor table | TCR distances and optional clusters |
| ImmGLIPH | `ImmunData` with TRB CDR3 amino-acid sequences | Groups of similar TRB receptors |
| IgBLAST | `ImmunData`, FASTA, or FASTA with an annotation table | AIRR table with germline-gene calls |
| RAxML-NG | `ImmunData` with aligned BCR sequences and family names | One phylogenetic tree per family |
| VDJdb | `ImmunData` passed to `data$annotate()` | Matching antigen-specificity records |

## tcrdist3

tcrdist3 calculates distances between T-cell receptors. It can also group
similar receptors.

Input:

- Pass an `ImmunData` dataset or one CSV, TSV, or TXT file.
- Each receptor needs a CDR3 amino-acid sequence, V gene, and J gene for every
  selected chain.
- For paired receptors in a text file, use one row per receptor and separate
  CDR3, V-gene, and J-gene columns for each chain.

Output:

- `results$distances` is the path to the NumPy distance matrix.
- `results$mapping` shows which identical receptors were calculated once.
- In cluster mode, `results$clusters` contains one row per receptor with
  `tcrdist_cluster_id` and `tcrdist_cluster_size`.
- `results$edges` contains receptor pairs whose distance is within `radius`.
- For `ImmunData` input, `results$annotate(idata)` returns a new dataset with
  the cluster columns.

```r
tool <- tool_tcrdist3(conf_tcrdist3(
  mode = "cluster",
  cpus = 4L,
  radius = 24L
))

run <- tool$run(input, output_dir = "analysis/tcrdist3")
results <- run$results()
clusters <- dplyr::collect(results$clusters)
```

## ImmGLIPH

ImmGLIPH groups TRB receptors that share sequence patterns.

Input:

- Pass an `ImmunData` dataset containing TRB receptors.
- By default, `junction_aa` contains the TRB CDR3 amino-acid sequence.
- V gene, donor, HLA, and clone-count columns are optional. Select them with
  `v_gene_col`, `patient_col`, `hla_col`, and `counts_col`.

Output:

- `results$groups` contains one row per grouped receptor and its
  `immgliph_group_id`.
- `results$connections` contains sequence pairs connected by ImmGLIPH.
- `results$cluster_properties` contains one row per group with its properties
  and scores.
- `results$mapping` contains the receptor columns sent to ImmGLIPH.
- `results$annotate(idata)` returns a new `ImmunData` dataset with the group
  column. Receptors without a group have `NA`.

```r
tool <- tool_immgliph(conf_immgliph(
  v_gene_col = "v_call",
  patient_col = "donor",
  seed = 1L
))

run <- tool$run(idata, output_dir = "analysis/immgliph")
results <- run$results()
groups <- dplyr::collect(results$groups)
```

## IgBLAST

IgBLAST assigns V, D, J, and optional constant-region genes to immune-receptor
nucleotide sequences. The germline database selected in `igblastr` determines
the organism and receptor type.

Choose one input:

1. `ImmunData`

   Pass an `ImmunData` dataset. `sequence_col` gives the column containing
   full-length, ungapped nucleotide sequences.

2. FASTA

   Pass one FASTA or gzip-compressed FASTA file. Every sequence must have a
   unique name and contain ungapped IUPAC nucleotide letters.

3. FASTA with an annotation table

   Pass two named file paths:

   ```r
   c(
     fasta = "filtered_contig.fasta",
     annotations = "filtered_contig_annotations.csv"
   )
   ```

   The annotation table may be CSV, TSV, or TXT. It must contain a unique
   `contig_id` or `sequence_id` column. Values in this column must match the
   FASTA sequence names.

Output:

- `results$airr` is the complete AIRR result table.
- `results$mapping` shows which input chain or FASTA sequence produced each
  AIRR `sequence_id`.
- For `ImmunData` input, `results$annotate(idata)` returns a new dataset with
  selected AIRR columns.
- For FASTA with an annotation table, `results$final` contains the original
  columns joined with the AIRR columns. `results$final_output` is the path to
  the same table as a TSV file.

```r
tool <- tool_igblast(conf_igblast(sequence_col = "sequence"))

run <- tool$run(idata, output_dir = "analysis/igblast")
results <- run$results()
airr <- dplyr::collect(results$airr)

aligned_idata <- results$annotate(idata)
```

In a FASTA header such as `>cell_1 description`, the sequence name is `cell_1`.
IgBLAST writes this name to the AIRR `sequence_id` column.

## RAxML-NG

RAxML-NG builds one phylogenetic tree per B-cell clone or family. It does not
identify clones, reconstruct germline sequences, or align sequences.

Input:

- Pass an `ImmunData` dataset.
- `family_col` gives the clone or family column. A separate tree is built for
  each value.
- `tip_col` gives the sequence name shown at each tree tip. Names must be
  unique within a family.
- `alignment_col` gives the aligned nucleotide sequence. All sequences in one
  family must have the same length.

Output:

- `results$families` shows whether each family analysis succeeded.
- `results$tree_index` gives the path to each successful tree file.
- `results$tip_map` contains the original tip names and the shorter names used
  by RAxML-NG.
- `results$read_trees()` returns one Newick tree per family.

```r
tool <- tool_raxml_ng(conf_raxml_ng(
  alignment_col = "clone_sequence_alignment",
  family_col = "clone_id"
))

run <- tool$run(idata, output_dir = "analysis/raxml-ng")
results <- run$results()
trees <- results$read_trees()
```

## VDJdb

VDJdb contains T-cell receptors with known antigen specificities.

Input:

- `data$load()` takes no biological input and reads the downloaded VDJdb table.
- `data$annotate(idata)` takes an `ImmunData` dataset.
- By default, matching uses the `junction_aa` and `locus` columns.

Output:

- `data$load()` returns the VDJdb table.
- `data$annotate(idata)` returns a new `ImmunData` dataset. Added columns start
  with `vdjdb_`.

```r
data <- data_vdjdb()
data$download(execute = TRUE)

vdjdb <- data$load()
annotated_idata <- data$annotate(idata)
```

## Files and data

By default, downloaded databases and tool environments are stored under
`~/immundapter_data`. Set `IMMUNDAPTER_DATA` to use another directory.

Each analysis has its own output directory containing the prepared input,
settings, logs, and result files.

## About

`immundapter` is part of the
[ImmunoMind ecosystem](https://github.com/immunomind) for immune-repertoire
data ingestion, storage, analysis, and annotation.

Author and maintainer: [Vadim I. Nazarov](https://orcid.org/0000-0003-3659-2709),
ImmunoMind.

For questions and support, use
[GitHub Issues](https://github.com/immunomind/immundapter/issues) or email
[`support@immunomind.com`](mailto:support@immunomind.com).
