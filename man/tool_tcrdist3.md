

# Compare and group T-cell receptors with tcrdist3

[**Source code**](https://github.com/immunomind/immundapter/tree/main/R/tcrdist3.R#L243)

## Description

tcrdist3 measures how similar T-cell receptors are from their CDR3
amino-acid sequences and V/J genes. A smaller distance means that two
receptors are more similar. The tool can return the distances alone or
use them to place similar receptors into groups.

The usual workflow is:

<pre>tool &lt;- tool_tcrdist3()
run &lt;- tool\$run(idata)
results &lt;- run\$results()
clusters &lt;- dplyr::collect(results\$clusters)
</pre>

## Usage

<pre><code class='language-R'>conf_tcrdist3(
  mode = c("cluster", "distance"),
  organism = "human",
  chains = "beta",
  cdr3_cols = c(beta = "junction_aa"),
  v_gene_cols = c(beta = "v_call"),
  j_gene_cols = c(beta = "j_call"),
  db_file = "alphabeta_gammadelta_db.tsv",
  store_all_cdr = TRUE,
  deduplicate = TRUE,
  cpus = 1L,
  radius = 24L,
  python = NULL,
  run_args = list()
)

tool_tcrdist3(config = conf_tcrdist3())
</code></pre>

## Arguments

<table role="presentation">
<tr>
<td style="white-space: nowrap; font-family: monospace; vertical-align: top">
<code id="mode">mode</code>
</td>
<td>
<code>“cluster”</code> groups similar receptors and returns their
distances. <code>“distance”</code> returns only the distance matrix and
its input mapping.
</td>
</tr>
<tr>
<td style="white-space: nowrap; font-family: monospace; vertical-align: top">
<code id="organism">organism</code>
</td>
<td>
Species used to interpret V/J genes: <code>“human”</code> or
<code>“mouse”</code>.
</td>
</tr>
<tr>
<td style="white-space: nowrap; font-family: monospace; vertical-align: top">
<code id="chains">chains</code>
</td>
<td>
One or two chains selected from <code>“alpha”</code>,
<code>“beta”</code>, <code>“gamma”</code>, and <code>“delta”</code>.
With two chains, their distances are added.
</td>
</tr>
<tr>
<td style="white-space: nowrap; font-family: monospace; vertical-align: top">
<code id="cdr3_cols">cdr3_cols</code>
</td>
<td>
Named vector giving the CDR3 amino-acid column for each selected chain,
for example <code>c(beta = “junction_aa”)</code>.
</td>
</tr>
<tr>
<td style="white-space: nowrap; font-family: monospace; vertical-align: top">
<code id="v_gene_cols">v_gene_cols</code>
</td>
<td>
Named vector giving the V-gene column for each selected chain. Gene
names must be recognized by the tcrdist3 database.
</td>
</tr>
<tr>
<td style="white-space: nowrap; font-family: monospace; vertical-align: top">
<code id="j_gene_cols">j_gene_cols</code>
</td>
<td>
Named vector giving the J-gene column for each selected chain.
</td>
</tr>
<tr>
<td style="white-space: nowrap; font-family: monospace; vertical-align: top">
<code id="db_file">db_file</code>
</td>
<td>
Germline reference supplied to tcrdist3. The default reference is
included with tcrdist3.
</td>
</tr>
<tr>
<td style="white-space: nowrap; font-family: monospace; vertical-align: top">
<code id="store_all_cdr">store_all_cdr</code>
</td>
<td>
Whether to save separate distance matrices for individual receptor
regions in addition to the combined matrices.
</td>
</tr>
<tr>
<td style="white-space: nowrap; font-family: monospace; vertical-align: top">
<code id="deduplicate">deduplicate</code>
</td>
<td>
Whether to calculate identical receptors once. Results are still
reported for every original receptor.
</td>
</tr>
<tr>
<td style="white-space: nowrap; font-family: monospace; vertical-align: top">
<code id="cpus">cpus</code>
</td>
<td>
Number of CPU cores used by tcrdist3.
</td>
</tr>
<tr>
<td style="white-space: nowrap; font-family: monospace; vertical-align: top">
<code id="radius">radius</code>
</td>
<td>
Largest distance that directly connects two receptors during clustering.
It is used only when <code>mode = “cluster”</code>.
</td>
</tr>
<tr>
<td style="white-space: nowrap; font-family: monospace; vertical-align: top">
<code id="python">python</code>
</td>
<td>
Python executable name or path. <code>NULL</code> selects the managed
tcrdist3 environment described above.
</td>
</tr>
<tr>
<td style="white-space: nowrap; font-family: monospace; vertical-align: top">
<code id="run_args">run_args</code>
</td>
<td>
Additional named arguments passed to <code>TCRrep</code>. This is for
advanced settings. Values must be valid JSON values. Arguments
controlled by immundapter cannot be changed here.
</td>
</tr>
<tr>
<td style="white-space: nowrap; font-family: monospace; vertical-align: top">
<code id="config">config</code>
</td>
<td>
Settings created by <code>conf_tcrdist3()</code>.
</td>
</tr>
</table>

## Details

tcrdist3 compares receptor sequences; it does not prove that receptors
bind the same antigen. Treat the distance and cluster threshold as
analysis choices that should be checked for the biological question and
dataset.

Every selected chain needs a CDR3 amino-acid sequence, a V-gene call,
and a J-gene call. The defaults expect human TRB data stored in
<code>junction_aa</code>, <code>v_call</code>, and <code>j_call</code>.

## Value

<code>conf_tcrdist3()</code> returns read-only tcrdist3 settings.
<code>tool_tcrdist3()</code> returns a reusable tcrdist3 tool.

## Choose distances or clusters

Use <code>mode = “cluster”</code>, the default, to calculate distances
and group receptors. Two receptors are connected when their distance is
at most <code>radius</code>. All receptors connected directly or through
other receptors become one cluster. Consequently, two receptors in the
same cluster can have a direct distance greater than
<code>radius</code>.

Use <code>mode = “distance”</code> when only the pairwise distance
matrix is needed. Distance-only results do not contain clusters and
cannot be added to <code>ImmunData</code> with
<code style="white-space: pre;">$annotate()</code>.

## Choose the receptor chains

<code>chains</code> selects the chains used in the distance:

<ul>
<li>

<code>chains = “beta”</code>, the default, compares TRB only. It can
select the beta chain from paired TRA/TRB <code>ImmunData</code>.

</li>
<li>

<code>chains = c(“alpha”, “beta”)</code> compares both sides of a paired
receptor. tcrdist3 calculates one distance per chain and adds them.

</li>
</ul>

Alpha, beta, gamma, and delta chains are supported.
<code>organism</code> must match the gene names in the data. Because
paired-chain distances are sums, a <code>radius</code> chosen for one
chain should not be assumed to work for two chains.

For paired <code>ImmunData</code>, the receptor schema tells the tool
which row belongs to each chain, so both chains can use the same feature
names:

<pre>conf_tcrdist3(
  chains = c("alpha", "beta"),
  cdr3_cols = c(alpha = "junction_aa", beta = "junction_aa"),
  v_gene_cols = c(alpha = "v_call", beta = "v_call"),
  j_gene_cols = c(alpha = "j_call", beta = "j_call")
)
</pre>

## Choose the input

<code>tool$run()</code> accepts:

<ul>
<li>

<strong>An <code>ImmunData</code> dataset.</strong> The selected columns
must be receptor features. The tool prepares one row per biological
receptor. Cluster results can be added back with
<code>results$annotate(idata)</code>.

</li>
<li>

<strong>A CSV, TSV, or TXT file.</strong> Each row must represent one
receptor. Paired receptors need separate CDR3, V-gene, and J-gene
columns for each chain.

</li>
</ul>

Rows missing a required value for a selected chain are excluded. For
text input, <code>imd_receptor_id</code> is used as the identifier when
present, followed by <code>cell_id</code>; otherwise, row numbers are
used.

## Run the tool

<code>tool_tcrdist3()</code> creates a reusable tool. Its main methods
are:

<ul>
<li>

<code>tool$help()</code>: print commands for installing the Python
environment;

</li>
<li>

<code>tool$check()</code>: check that Python and tcrdist3 are ready;

</li>
<li>

<code>tool$run(input, config = NULL, output_dir = NULL)</code>: run one
analysis.

</li>
</ul>

The tool keeps its default settings in <code>tool$config</code>.
Supplying <code>config</code> to
<code style="white-space: pre;">$run()</code> replaces those settings
for that run only. The completed run stores its output directory in
<code style="white-space: pre;">$path</code> and its saved files in
<code style="white-space: pre;">$artifacts</code>.

## Read the results

Call <code>run$results()</code> after a successful run. Both modes
provide:

<ul>
<li>

<code>results$distances</code>: the path to the complete NumPy distance
matrix;

</li>
<li>

<code>results$mapping</code>: the link between input receptors and rows
of that matrix.

</li>
</ul>

Cluster mode also provides:

<ul>
<li>

<code>results$clusters</code>: one cluster number and cluster size per
receptor;

</li>
<li>

<code>results$edges</code>: receptor pairs whose distance is at most
<code>radius</code>;

</li>
<li>

<code>results$annotate(idata)</code>: a new <code>ImmunData</code>
dataset containing <code>tcrdist_cluster_id</code> and
<code>tcrdist_cluster_size</code>.

</li>
</ul>

Result tables stay on disk and are evaluated with duckplyr. Use
<code>dplyr::collect()</code> when a table is small enough to bring into
R memory. The distance matrix is a NumPy <code>.npy</code> file, not an
R table.

## Calculation details

With <code>deduplicate = TRUE</code>, identical selected receptors are
calculated once and then matched back to every original receptor.
<code>results$mapping</code> records this relationship. Set
<code>store_all_cdr = FALSE</code> when separate region-specific
distance matrices are not needed.

## Python setup

By default, immundapter uses a dedicated Python environment at
<code>environments/tcrdist3</code> under
<code>tools::R_user_dir(“immundapter”, “data”)</code>. Call
<code>tool$help()</code> to print the installation commands and
<code>tool$check()</code> to verify the environment.
<code style="white-space: pre;">$help()</code> prints instructions; it
does not install software.

Set <code>IMMUNDAPTER_DATA</code> before creating the configuration to
use another data directory, or supply <code>python</code> to use a
specific Python executable.

## See Also

immundata::ImmunData,
<a href="https://tcrdist3.readthedocs.io/en/latest/inputs.html">tcrdist3
input</a>, <a href="https://tcrdist3.readthedocs.io/en/latest/">tcrdist3
documentation</a>

## Examples

``` r
library("immundapter")


# Group similar TRB receptors in the bundled example.
config <- conf_tcrdist3(
  store_all_cdr = FALSE,
  cpus = 1L,
  radius = 24L
)
tool <- tool_tcrdist3(config)
stopifnot(tool$check())

input_file <- system.file(
  "extdata", "tcrdist3", "tcr-cluster-simple.tsv",
  package = "immundapter"
)
run <- tool$run(
  input_file,
  output_dir = tempfile("tcrdist3-file-")
)
results <- run$results()
clusters <- dplyr::collect(results$clusters)

# Add the cluster assignments to ImmunData.
snapshot <- system.file(
  "extdata", "tcrdist3", "immundata-tcr-cluster-simple",
  package = "immundapter"
)
idata <- immundata::read_immundata(snapshot, verbose = FALSE)
idata_run <- tool$run(
  idata,
  output_dir = tempfile("tcrdist3-immundata-")
)
annotated_idata <- idata_run$results()$annotate(idata)

# Calculate distances without creating clusters.
distance_tool <- tool_tcrdist3(conf_tcrdist3(
  mode = "distance", store_all_cdr = FALSE
))
distance_run <- distance_tool$run(
  input_file,
  output_dir = tempfile("tcrdist3-distance-")
)
distance_file <- distance_run$results()$distances
```
