

# VDJdb data source

[**Source code**](https://github.com/immunomind/immundapter/tree/main/R/vdjdb.R#L62)

## Description

VDJdb is a reference database of T-cell receptors with known antigen
specificities. <code>data_vdjdb()</code> manages a local copy of the
slim VDJdb table and can add matching VDJdb records to
immundata::ImmunData.

## Usage

<pre><code class='language-R'>data_vdjdb(output_dir = immundapter_data_path("data"))
</code></pre>

## Arguments

<table role="presentation">
<tr>
<td style="white-space: nowrap; font-family: monospace; vertical-align: top">
<code id="output_dir">output_dir</code>
</td>
<td>
Directory used to store VDJdb. Defaults to the <code>data</code>
subdirectory of <code>IMMUNDAPTER_DATA</code> when set, or of
<code>tools::R_user_dir(“immundapter”, “cache”)</code> otherwise.
</td>
</tr>
</table>

## Value

A VDJdb data source.

## Input and output

<code>data \<- data_vdjdb()</code> gives access to these commands:

<ul>
<li>

<strong><code>data$download()</code>:</strong> Print the download
command and destination without downloading anything.

</li>
<li>

<strong><code>data$download(execute = TRUE)</code>:</strong> Download
the latest VDJdb release. Use <code>force = TRUE</code> to replace an
existing copy.

</li>
<li>

<strong><code>data$check()</code>:</strong> Report whether VDJdb has
been downloaded.

</li>
<li>

<strong><code>data$load()</code>:</strong> Return the VDJdb table. The
table stays on disk until you use <code>dplyr::collect()</code>. Use
<code>data$load(output = “df”)</code> to return a base R data frame
instead.

</li>
<li>

<strong><code>data$annotate(idata)</code>:</strong> Match receptors in
an <code>ImmunData</code> dataset to VDJdb and return a new dataset with
VDJdb columns. By default, matching uses CDR3 amino-acid sequence and
receptor locus.

</li>
<li>

<strong><code>data$citation()</code>:</strong> Return the VDJdb
publication as an R citation.

</li>
</ul>

For the default matching, <code>ImmunData</code> must contain
<code>junction_aa</code> and <code>locus</code>. These are matched to
the VDJdb columns <code>cdr3</code> and <code>gene</code>:
<code>c(junction_aa = “cdr3”, locus = “gene”)</code>. Use the
<code>by</code> argument to choose different matching columns. Added
VDJdb columns start with <code>vdjdb\_</code>, which makes them easy to
distinguish from the original columns.

## See Also

<a href="https://github.com/antigenomics/vdjdb-db">VDJdb database</a>

## Examples

``` r
library("immundapter")


data <- data_vdjdb(output_dir = tempfile("vdjdb-cache-"))

# First call: print instructions. Second call: perform the download.
data$download()
data$download(execute = TRUE)
data$check()

vdjdb_table <- data$load()
dplyr::collect(vdjdb_table)

input_file <- system.file(
  "extdata", "tcrdist3", "tcr-cluster-simple.tsv",
  package = "immundapter"
)
idata <- immundata::read_repertoires(
  path = input_file,
  schema = immundata::make_receptor_schema(
    features = c("v_call", "j_call", "junction_aa"), chains = "TRB"
  ),
  barcode_col = "cell_id", locus_col = "locus",
  umi_col = "duplicate_count",
  output_folder = tempfile("vdjdb-immundata-")
)
annotated_idata <- data$annotate(idata)
```
