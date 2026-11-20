

# Assign germline genes to raw receptor sequences with IgBLAST

[**Source code**](https://github.com/immunomind/immundapter/tree/main/R/igblast.R#L212)

## Description

IgBLAST identifies the germline genes that most likely produced each
immune receptor sequence. It produces V, D, J, and, when available, C
gene calls in an AIRR table that can be imported with
<code>immundata::read_repertoires()</code>.

The usual workflow is:

``` r
run <- tool_igblast()$run(c(
  fasta = "sequences.fasta",
  annotations = "contig_annotations.csv"
))
idata <- immundata::read_repertoires(
  path = run$results()$final_output,
  schema = immundata::make_receptor_schema(
    features = c("v_call", "j_call", "junction_aa"),
    chains = c("IGH", "IGK|IGL")
  ),
  barcode_col = "barcode",
  locus_col = "locus",
  umi_col = "umis",
  preprocess = immundata::make_default_preprocessing("10x"),
  rename_columns = NULL
)
```

## Usage

<pre><code class='language-R'>conf_igblast(num_threads = NULL, run_args = list())

tool_igblast(config = conf_igblast())
</code></pre>

## Arguments

<table role="presentation">
<tr>
<td style="white-space: nowrap; font-family: monospace; vertical-align: top">
<code id="num_threads">num_threads</code>
</td>
<td>
Number of computer threads available to IgBLAST. <code>NULL</code> uses
the IgBLAST default.
</td>
</tr>
<tr>
<td style="white-space: nowrap; font-family: monospace; vertical-align: top">
<code id="run_args">run_args</code>
</td>
<td>
Additional named arguments passed through <code>…</code> in
<code>igblastr::igblastn()</code>. This is intended for advanced
settings. Arguments required by the immundapter workflow cannot be
changed here.
</td>
</tr>
<tr>
<td style="white-space: nowrap; font-family: monospace; vertical-align: top">
<code id="config">config</code>
</td>
<td>
Settings created by <code>conf_igblast()</code>.
</td>
</tr>
</table>

## Details

IgBLAST compares observed receptor sequences with a germline database.
The database selected through <code>igblastr</code> determines the
species and receptor type. For example, select a human BCR database
before analysing human immunoglobulin sequences.

Query sequences must be <strong>nucleotide sequences without alignment
gaps</strong>. Supply the complete assembled V(D)J sequence when
possible. Amino-acid sequences and sequences that already contain IMGT
dots are not valid input. A CDR3 nucleotide sequence alone usually does
not contain enough information for useful V- and J-gene annotation.

## Value

<code>conf_igblast()</code> returns read-only IgBLAST settings.
<code>tool_igblast()</code> returns a reusable IgBLAST tool.

## Before the first analysis

Install IgBLAST and one germline database with <code>igblastr</code>.
Select the database in the same R session in which you run the analysis:

``` r
igblastr::install_igblast()
igblastr::install_IMGT_germline_db(
  "<release>",
  organism = "Homo sapiens",
  tcr.db = FALSE
)
igblastr::use_germline_db("<db-name>")
```

Use <code>tool$help()</code> for setup commands and
<code>tool$check()</code> to confirm that the selected database and
IgBLAST installation are ready.

## Run the tool

<code>tool_igblast()</code> creates a reusable tool. Its main methods
are:

<ul>
<li>

<code>tool$help()</code>: print installation and database setup
instructions;

</li>
<li>

<code>tool$check()</code>: check that IgBLAST and a germline database
are ready;

</li>
<li>

<code>tool$run(input, config = NULL, output_dir = NULL)</code>: run one
analysis.

</li>
</ul>

The tool keeps its default settings in <code>tool$config</code>.
Supplying <code>config</code> to
<code style="white-space: pre;">$run()</code> replaces those settings
for that run only.

## Choose the input

<code>tool$run()</code> accepts two forms of input:

<ul>
<li>

<strong>A FASTA file.</strong> Pass a path ending in <code>.fasta</code>
or <code>.fasta.gz</code>. Use this when no annotation table needs to be
preserved.

</li>
<li>

<strong>A FASTA file and its annotation table.</strong> Pass
<code>c(fasta = “sequences.fasta”, annotations =
“annotations.csv”)</code>. Use this for data such as 10x contigs when
barcodes, UMI counts, and sample columns must remain beside the new
IgBLAST annotations.

</li>
</ul>

## Read the results

<code>tool$run()</code> saves the analysis and returns a completed run
object. Its <code style="white-space: pre;">$path</code> is the output
directory and <code style="white-space: pre;">$artifacts</code> lists
the saved files. Call <code>run$results()</code> to access the
biological results. Every run provides:

<ul>
<li>

<code>results$airr</code>: the complete AIRR annotation table;

</li>
<li>

<code>results$mapping</code>: the link between input records and AIRR
<code>sequence_id</code>;

</li>
<li>

<code>results$raw_output</code>: the path to the saved AIRR TSV file.

</li>
</ul>

The tables are disk-backed. Use <code>dplyr::collect()</code> only when
you want to bring a table into R memory.

With paired FASTA and annotation input, <code>results$final</code>
contains the source columns joined with all AIRR columns.
<code>results$final_output</code> is the path to the same table as a TSV
file. Pass this file to <code>immundata::read_repertoires()</code> to
create an <code>ImmunData</code> dataset. Choose preprocessing that
retains <code>sequence_alignment</code> when subsequent analyses need
IMGT-gapped alignments and coordinates.

## IMGT gaps and coordinates

<code>igblastr</code> builds the searchable IgBLAST database from
<strong>ungapped</strong> germline sequences. It also keeps the original
<strong>IMGT-gapped V references</strong>. These two forms have
different V-gene coordinates.

For a selected IMGT database, immundapter automatically uses both forms
from the <strong>same installed database release</strong>. IgBLAST first
performs the search against the ungapped database. immundapter then uses
the matching gapped V allele to place IMGT dots in
<code>sequence_alignment</code> and convert
<code>v_germline_start</code> and <code>v_germline_end</code> to the
IMGT-gapped coordinate system. This follows the normalization used by
Change-O for germline reconstruction.

Do not download a newer gapped FASTA separately: allele content and gap
positions must match the database used by IgBLAST. If the selected IMGT
database does not contain its original gapped V FASTA files, reinstall
that database with <code>igblastr::install_IMGT_germline_db()</code>.

Normalization uses the first allele in <code>v_call</code>. If the
allele or the exact aligned V end cannot be matched to the selected
reference, that row keeps the original IgBLAST alignment and
coordinates. Non-IMGT databases are also returned in their original
IgBLAST coordinate system. D- and J-region coordinates are not converted
to IMGT-gapped V coordinates.

## Input details

FASTA names must be unique. In <code style="white-space: pre;">\>cell_1
description</code>, the name is <code>cell_1</code>. For paired input,
the annotation table must contain a unique <code>contig_id</code> or
<code>sequence_id</code> column whose values match the FASTA names. CSV,
TSV, TXT, and gzip-compressed annotation files are accepted.

## See Also

<a href="https://bioconductor.org/packages/release/bioc/vignettes/igblastr/inst/doc/igblastr_overview.html">igblastr
overview</a>

## Examples

``` r
library("immundapter")


# Select an installed human BCR IMGT database.
databases <- igblastr::list_germline_dbs(names.only = TRUE)
databases <- grep(
  "Homo_sapiens.IGH+IGK+IGL",
  databases,
  value = TRUE,
  fixed = TRUE
)
databases <- grep("^IMGT-", databases, value = TRUE)
stopifnot(length(databases) > 0L)
suppressMessages(
  igblastr::use_germline_db(sort(databases)[[length(databases)]])
)

config <- conf_igblast()
tool <- tool_igblast(config)
stopifnot(tool$check())

antitnf_dir <- system.file("extdata", "antitnf", package = "immundapter")
antitnf_fasta <- file.path(antitnf_dir, "bcr_contigs.fasta.gz")
antitnf_annotations <- file.path(
  antitnf_dir,
  "contig_annotations.csv.gz"
)

# Analyse a FASTA file directly.
fasta_run <- tool$run(
  antitnf_fasta,
  output_dir = tempfile("igblast-fasta-")
)
fasta_airr <- dplyr::collect(fasta_run$results()$airr)

# Preserve the matching 10x annotation table and import the AIRR result.
paired_run <- tool$run(
  c(
    fasta = antitnf_fasta,
    annotations = antitnf_annotations
  ),
  output_dir = tempfile("igblast-10x-")
)
idata <- immundata::read_repertoires(
  path = paired_run$results()$final_output,
  schema = immundata::make_receptor_schema(
    features = c("v_call", "j_call", "junction_aa"),
    chains = c("IGH", "IGK|IGL")
  ),
  barcode_col = "barcode",
  locus_col = "locus",
  umi_col = "umis",
  preprocess = immundata::make_default_preprocessing("10x"),
  rename_columns = NULL,
  output_folder = tempfile("antitnf-immundata-")
)
```
