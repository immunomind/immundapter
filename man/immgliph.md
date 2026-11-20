

# Group similar T-cell receptors with ImmGLIPH

## Description

ImmGLIPH groups T-cell receptor beta chains that share sequence patterns
and may recognize related antigens. The analysis uses CDR3 amino-acid
sequences and can also use V genes, donor names, HLA types, and clone
counts.

<code>ImmGliphToolConfig</code> stores the validated settings for an
ImmGLIPH analysis. It contains no repertoire data and no results. Create
it with <code>conf_immgliph()</code> and pass it to
<code>tool_immgliph()</code>.

## Usage

<pre><code class='language-R'>conf_immgliph(
  method = "gliph2",
  reference = "human_v2.0_CD48",
  cdr3_col = "junction_aa",
  v_gene_col = NULL,
  patient_col = NULL,
  hla_col = NULL,
  counts_col = NULL,
  n_cores = 1L,
  seed = NULL,
  run_args = list()
)

tool_immgliph(config = conf_immgliph())
</code></pre>

## Arguments

<table role="presentation">
<tr>
<td style="white-space: nowrap; font-family: monospace; vertical-align: top">
<code id="method">method</code>
</td>
<td>
ImmGLIPH analysis preset: <code>“gliph2”</code>, <code>“gliph1”</code>,
or <code>“custom”</code>.
</td>
</tr>
<tr>
<td style="white-space: nowrap; font-family: monospace; vertical-align: top">
<code id="reference">reference</code>
</td>
<td>
Name of the reference repertoire included with ImmGLIPH.
</td>
</tr>
<tr>
<td style="white-space: nowrap; font-family: monospace; vertical-align: top">
<code id="cdr3_col">cdr3_col</code>
</td>
<td>
Name of the column containing TRB CDR3 amino-acid sequences.
</td>
</tr>
<tr>
<td style="white-space: nowrap; font-family: monospace; vertical-align: top">
<code id="v_gene_col">v_gene_col</code>
</td>
<td>
Optional name of the TRB V-gene column.
</td>
</tr>
<tr>
<td style="white-space: nowrap; font-family: monospace; vertical-align: top">
<code id="patient_col">patient_col</code>
</td>
<td>
Optional name of the donor or sample column.
</td>
</tr>
<tr>
<td style="white-space: nowrap; font-family: monospace; vertical-align: top">
<code id="hla_col">hla_col</code>
</td>
<td>
Optional name of the HLA column.
</td>
</tr>
<tr>
<td style="white-space: nowrap; font-family: monospace; vertical-align: top">
<code id="counts_col">counts_col</code>
</td>
<td>
Optional name of the clone-count column.
</td>
</tr>
<tr>
<td style="white-space: nowrap; font-family: monospace; vertical-align: top">
<code id="n_cores">n_cores</code>
</td>
<td>
Positive integer number of CPU cores used by ImmGLIPH.
</td>
</tr>
<tr>
<td style="white-space: nowrap; font-family: monospace; vertical-align: top">
<code id="seed">seed</code>
</td>
<td>
Optional integer random seed for a reproducible ImmGLIPH run.
</td>
</tr>
<tr>
<td style="white-space: nowrap; font-family: monospace; vertical-align: top">
<code id="run_args">run_args</code>
</td>
<td>
Additional named arguments passed to <code>immGLIPH::runGLIPH()</code>.
Use this only for advanced ImmGLIPH settings. Arguments controlled by
immundapter cannot be changed here.
</td>
</tr>
<tr>
<td style="white-space: nowrap; font-family: monospace; vertical-align: top">
<code id="config">config</code>
</td>
<td>
An ImmGliphToolConfig created by <code>conf_immgliph()</code>.
</td>
</tr>
</table>

## Value

<code>conf_immgliph()</code> returns a read-only ImmGliphToolConfig.
<code>tool_immgliph()</code> returns a reusable ImmGliphTool.

## Input to <code>tool$run()</code>

<code>tool$run()</code> accepts one immundata::ImmunData dataset. It
does not accept a text file directly. The dataset must contain TRB
receptors. For paired TRA/TRB data, only the TRB chain is used.

<code>cdr3_col</code> must name a column containing TRB CDR3 amino-acid
sequences. Missing and empty CDR3 values are skipped. The other input
columns are optional: set <code>v_gene_col</code>,
<code>patient_col</code>, <code>hla_col</code>, or
<code>counts_col</code> only when the analysis should use those data.
Every selected column must exist in the <code>ImmunData</code> dataset.

## Output from <code>tool$run()</code>

<code>tool$run()</code> returns an ImmGliphRun. It records where the
results were saved. To read the result tables, first run <code>results
\<- run$results()</code>.

ImmGliphResults provides:

<ul>
<li>

<strong><code>results$groups</code>:</strong> One row for each receptor
assigned to a group. The <code>immgliph_group_id</code> column gives the
group. Receptors without a group are absent. Use
<code>dplyr::collect()</code> to read the table into R.

</li>
<li>

<strong><code>results$connections</code>:</strong> Sequence pairs
connected by ImmGLIPH and the reason for each connection.

</li>
<li>

<strong><code>results$cluster_properties</code>:</strong> One row per
group with the properties and scores calculated by ImmGLIPH.

</li>
<li>

<strong><code>results$mapping</code>:</strong> The receptor columns that
were passed to ImmGLIPH.

</li>
<li>

<strong><code>results$raw_output</code>:</strong> Path to the original
ImmGLIPH <code>.rds</code> file.

</li>
<li>

<strong><code>results$annotate(idata)</code>:</strong> A new
<code>ImmunData</code> dataset with an <code>immgliph_group_id</code>
receptor column. Receptors without a group have <code>NA</code>.

</li>
</ul>

## Super class

<code>AbstractConfig</code> -\> <code>ImmGliphToolConfig</code>

## Public fields

<dl>
<dt>
<code>method</code>
</dt>
<dd>
ImmGLIPH preset: <code>“gliph2”</code>, <code>“gliph1”</code>, or
<code>“custom”</code>.
</dd>
<dt>
<code>reference</code>
</dt>
<dd>
Built-in ImmGLIPH reference repertoire name.
</dd>
<dt>
<code>cdr3_col</code>
</dt>
<dd>
Column containing CDR3 amino-acid sequences.
</dd>
<dt>
<code>v_gene_col</code>
</dt>
<dd>
Optional V-gene column.
</dd>
<dt>
<code>patient_col</code>
</dt>
<dd>
Optional donor or sample column.
</dd>
<dt>
<code>hla_col</code>
</dt>
<dd>
Optional HLA column.
</dd>
<dt>
<code>counts_col</code>
</dt>
<dd>
Optional clone-count column.
</dd>
<dt>
<code>n_cores</code>
</dt>
<dd>
Number of cores passed to ImmGLIPH.
</dd>
<dt>
<code>seed</code>
</dt>
<dd>
Optional random seed used around the ImmGLIPH call.
</dd>
<dt>
<code>run_args</code>
</dt>
<dd>
Additional named arguments for <code>immGLIPH::runGLIPH()</code>.
</dd>
</dl>

## Methods

<h4>
Public methods
</h4>
<ul>
<li>

<a href="#method-ImmGliphToolConfig-initialize"><code>ImmGliphToolConfig$new()</code></a>

</li>
</ul>
<details open>
<summary>
Inherited methods
</summary>
<ul>
<li>
<span class="pkg-link" data-pkg="immundapter"
data-topic="AbstractConfig"
data-id="as_list"><a href='../../immundapter/html/AbstractConfig.html#method-AbstractConfig-as_list'><code>AbstractConfig$as_list()</code></a></span>
</li>
</ul>
</details>
<hr>

<a id="method-ImmGliphToolConfig-initialize"></a>

<h4>
<code>ImmGliphToolConfig$new()</code>
</h4>

Create checked ImmGLIPH settings. Users should call
<code>conf_immgliph()</code> instead.

<h5>
Usage
</h5>

<pre>ImmGliphToolConfig\$new(values)</pre>

<h5>
Arguments
</h5>

<dl>
<dt>
<code>values</code>
</dt>
<dd>
A named list of ImmGLIPH settings.
</dd>
</dl>

## See Also

<a href="https://bioc.r-universe.dev/immGLIPH/doc/manual.html#runGLIPH">ImmGLIPH
manual</a>

## Examples

``` r
library("immundapter")


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
  output_folder = tempfile("immgliph-immundata-")
)
config <- conf_immgliph(
  v_gene_col = "v_call",
  n_cores = 4L,
  seed = 1L
)
tool <- tool_immgliph(config)
tool$help()

if (tool$check()) {
  run <- tool$run(idata, output_dir = tempfile("immgliph-run-"))
  results <- run$results()

  groups <- dplyr::collect(results$groups)
  group_scores <- dplyr::collect(results$cluster_properties)
  grouped_idata <- results$annotate(idata)
}
```
