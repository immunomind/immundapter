

# Infer family trees with RAxML-NG

## Description

RAxML-NG builds a phylogenetic tree for each sequence family. This tool
takes aligned nucleotide sequences from immundata::ImmunData and runs
each family independently. It does not define families, align sequences,
or reconstruct germlines. Complete those steps before running RAxML-NG.

<code>RaxmlNgToolConfig</code> stores validated settings for one type of
family-level RAxML-NG analysis. It contains no sequences and no trees.
Create it with <code>conf_raxml_ng()</code> and pass it to
<code>tool_raxml_ng()</code>. The settings cannot be changed after
creation; call <code>conf_raxml_ng()</code> again to use different
settings.

## Usage

<pre><code class='language-R'>conf_raxml_ng(
  alignment_col = "clone_sequence_alignment",
  family_col = "clone_id",
  tip_col = immundata::imd_schema("receptor"),
  locus = "IGH",
  analysis = c("search", "fast"),
  model = "GTR+G",
  seed = NULL,
  threads = NULL,
  ancestral = FALSE,
  partition_col = NULL,
  brlen = NULL
)

tool_raxml_ng(config = conf_raxml_ng())
</code></pre>

## Arguments

<table role="presentation">
<tr>
<td style="white-space: nowrap; font-family: monospace; vertical-align: top">
<code id="alignment_col">alignment_col</code>
</td>
<td>
Name of the column containing the aligned nucleotide sequence for each
tree tip. This is not IgBLAST’s AIRR <code>sequence_alignment</code>
field, which aligns each sequence to a germline rather than to the other
sequences in its family.
</td>
</tr>
<tr>
<td style="white-space: nowrap; font-family: monospace; vertical-align: top">
<code id="family_col">family_col</code>
</td>
<td>
Column containing B-cell clone, lineage, or family names or numbers.
Each value is analyzed independently.
</td>
</tr>
<tr>
<td style="white-space: nowrap; font-family: monospace; vertical-align: top">
<code id="tip_col">tip_col</code>
</td>
<td>
Name of the column containing the sequence name used at each tree tip.
The default is the <code>ImmunData</code> receptor ID column.
</td>
</tr>
<tr>
<td style="white-space: nowrap; font-family: monospace; vertical-align: top">
<code id="locus">locus</code>
</td>
<td>
Receptor locus to use. It must agree with the receptor chains declared
by the <code>ImmunData</code> schema.
</td>
</tr>
<tr>
<td style="white-space: nowrap; font-family: monospace; vertical-align: top">
<code id="analysis">analysis</code>
</td>
<td>
Tree-inference mode. <code>“search”</code> performs the standard
maximum-likelihood search; <code>“fast”</code> performs the faster
search.
</td>
</tr>
<tr>
<td style="white-space: nowrap; font-family: monospace; vertical-align: top">
<code id="model">model</code>
</td>
<td>
Non-empty RAxML-NG nucleotide model string. The same model is used for
every partition in a partitioned analysis.
</td>
</tr>
<tr>
<td style="white-space: nowrap; font-family: monospace; vertical-align: top">
<code id="seed">seed</code>
</td>
<td>
Optional positive integer random seed. <code>NULL</code> retains the
RAxML-NG default.
</td>
</tr>
<tr>
<td style="white-space: nowrap; font-family: monospace; vertical-align: top">
<code id="threads">threads</code>
</td>
<td>
Optional positive integer thread count. <code>NULL</code> retains the
RAxML-NG automatic default.
</td>
</tr>
<tr>
<td style="white-space: nowrap; font-family: monospace; vertical-align: top">
<code id="ancestral">ancestral</code>
</td>
<td>
Whether to run RAxML-NG ancestral-state reconstruction on every
successfully inferred family tree.
</td>
</tr>
<tr>
<td style="white-space: nowrap; font-family: monospace; vertical-align: top">
<code id="partition_col">partition_col</code>
</td>
<td>
Optional column containing one partition specification per family.
</td>
</tr>
<tr>
<td style="white-space: nowrap; font-family: monospace; vertical-align: top">
<code id="brlen">brlen</code>
</td>
<td>
Optional branch-length linkage mode: <code>“linked”</code>,
<code>“scaled”</code>, or <code>“unlinked”</code>. This setting requires
<code>partition_col</code>; partitioned runs default to
<code>“scaled”</code>.
</td>
</tr>
<tr>
<td style="white-space: nowrap; font-family: monospace; vertical-align: top">
<code id="config">config</code>
</td>
<td>
A RaxmlNgToolConfig created by <code>conf_raxml_ng()</code>.
</td>
</tr>
</table>

## Value

<code>conf_raxml_ng()</code> returns a read-only RaxmlNgToolConfig.
<code>tool_raxml_ng()</code> returns a reusable RaxmlNgTool.

## Input to <code>tool$run()</code>

<code>tool$run()</code> accepts one immundata::ImmunData dataset. It
does not accept a FASTA file directly. The dataset must contain three
columns:

<ul>
<li>

<strong><code>family_col</code>:</strong> Clone, lineage, or family
name. A separate tree is built for each value in this column.

</li>
<li>

<strong><code>tip_col</code>:</strong> Sequence name shown at the tips,
or leaves, of the tree. Each sequence within a family must have a
different name.

</li>
<li>

<strong><code>alignment_col</code>:</strong> Aligned nucleotide
sequence. All sequences in one family must have the same length.

</li>
</ul>

All three values must be present and non-empty. Repeated rows with the
same family, tip, and aligned sequence are used once. A tip cannot have
two different sequences in the same family.

<code>additional_tips</code> optionally supplies backend-only tips that
are appended after the observed tips are extracted. It must be a data
frame or duckplyr table with <code>family_id</code>,
<code>tip_id</code>, and <code>aligned_sequence</code> columns. Its
family values must already occur in the extracted input. In a
partitioned analysis, an optional <code>partition_spec</code> column
must agree with the observed family specification; when omitted, that
specification is inherited. Additional tips are never added to or joined
into <code>ImmunData</code>.

<code>outgroups</code> optionally identifies one or more existing tips
per family. It must be a data frame or duckplyr table with
<code>family_id</code> and <code>tip_id</code> columns. The identified
tips may come from <code>ImmunData</code> or
<code>additional_tips</code>. RAxML-NG output is rooted by the adapter,
which writes a rooted copy of each successful tree while preserving the
raw RAxML-NG tree.

<code>locus</code> selects the receptor chain used for the tree. With
paired receptors, only the selected chain is used. The selected locus
must agree with the chains stored in the <code>ImmunData</code> dataset.

If <code>partition_col</code> is set, it must contain one consistent
partition specification per family in the form
<code>name=start-end</code>, with partitions separated by semicolons.
For example,
<code style="white-space: pre;">IGH=1-300;IGK=301-540</code> covers a
540-site alignment. Partitions must be unique, non-overlapping, and
cover every alignment site.

## Analysis modes

The tool supports two RAxML-NG tree-inference commands:

<ul>
<li>

<strong><code>“search”</code>:</strong> Standard maximum-likelihood tree
search.

</li>
<li>

<strong><code>“fast”</code>:</strong> Faster tree search.

</li>
</ul>

Set <code>ancestral = TRUE</code> to run RAxML-NG ancestral-state
reconstruction after each successful family tree. This produces a
node-labelled tree, reconstructed internal-node sequences, and marginal
state probabilities. A failed ancestral stage does not discard the
successful family tree. Other modes that use existing trees or produce
other types of output, including <code>evaluate</code>,
<code>bootstrap</code>, and <code>support</code>, are not part of this
tool. The default nucleotide model is <code>GTR+G</code>. Set
<code>seed</code> and <code>threads</code> for reproducible execution,
or omit them to retain RAxML-NG’s defaults. Partitioned analyses default
to scaled branch-length linkage.

## Output from <code>tool$run()</code>

<code>tool$run()</code> returns a RaxmlNgRun. It records where the
results were saved. To read the trees, first run <code>results \<-
run$results()</code>.

RaxmlNgResults provides:

<ul>
<li>

<strong><code>results$families</code>:</strong> One row per family,
showing whether the analysis succeeded and where its files were saved.
Use <code>dplyr::collect()</code> to read the table into R.

</li>
<li>

<strong><code>results$tree_index</code>:</strong> The families that
produced a tree and the path to each tree file.

</li>
<li>

<strong><code>results$failed_families</code>:</strong> Families that
failed, with their exit status and diagnostic message.

</li>
<li>

<strong><code>results$tip_map</code>:</strong> The original value from
<code>tip_col</code> and the shorter tip label used in the RAxML-NG
files.

</li>
<li>

<strong><code>results$outgroups</code>:</strong> The original and
internal outgroup tip names.

</li>
<li>

<strong><code>results$ancestral_index</code>:</strong> Paths to
successful ancestral-state files.

</li>
<li>

<strong><code>results$failed_ancestral</code>:</strong> Families whose
optional ancestral stage failed, with diagnostics.

</li>
<li>

<strong><code>results$read_trees()</code>:</strong> One Newick tree per
family, returned as a named character vector. When outgroups were
supplied, these are rooted copies.

</li>
<li>

<strong><code>results$read_ancestral_trees()</code>:</strong>
Node-labelled ancestral trees.

</li>
<li>

<strong><code>results$read_ancestral_states()</code>:</strong> A
duckplyr table of reconstructed internal-node sequences.

</li>
</ul>

Tree files remain in the run directory. Trees cannot be added as columns
to <code>ImmunData</code>. Use
<code style="white-space: pre;">$read_trees()</code>,
<code style="white-space: pre;">$tree_index</code>, or
<code style="white-space: pre;">$tip_map</code> for downstream
phylogenetic analysis.

Each family has its own directory containing the prepared FASTA file,
optional partition file, RAxML-NG output, and log files. A failed family
does not discard successful results. The returned run has status
<code>“partial”</code> when one or more tree or requested ancestral
stages fail.

## Super class

<code>AbstractConfig</code> -\> <code>RaxmlNgToolConfig</code>

## Public fields

<dl>
<dt>
<code>alignment_col</code>
</dt>
<dd>
Column containing family-aligned nucleotides.
</dd>
<dt>
<code>family_col</code>
</dt>
<dd>
Column containing the family name or number.
</dd>
<dt>
<code>tip_col</code>
</dt>
<dd>
Column containing the sequence name used at each tree tip.
</dd>
<dt>
<code>locus</code>
</dt>
<dd>
Receptor locus analyzed.
</dd>
<dt>
<code>analysis</code>
</dt>
<dd>
RAxML-NG inference mode.
</dd>
<dt>
<code>model</code>
</dt>
<dd>
RAxML-NG nucleotide model.
</dd>
<dt>
<code>seed</code>
</dt>
<dd>
Optional RAxML-NG random seed.
</dd>
<dt>
<code>threads</code>
</dt>
<dd>
Optional RAxML-NG thread count.
</dd>
<dt>
<code>ancestral</code>
</dt>
<dd>
Whether to reconstruct internal-node ancestral states.
</dd>
<dt>
<code>partition_col</code>
</dt>
<dd>
Optional partition-specification column.
</dd>
<dt>
<code>brlen</code>
</dt>
<dd>
Optional partition branch-length linkage mode.
</dd>
</dl>

## Methods

<h4>
Public methods
</h4>
<ul>
<li>

<a href="#method-RaxmlNgToolConfig-initialize"><code>RaxmlNgToolConfig$new()</code></a>

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

<a id="method-RaxmlNgToolConfig-initialize"></a>

<h4>
<code>RaxmlNgToolConfig$new()</code>
</h4>

Create checked RAxML-NG settings. Users should call
<code>conf_raxml_ng()</code> instead.

<h5>
Usage
</h5>

``` r
RaxmlNgToolConfig$new(values)
```

<h5>
Arguments
</h5>

<dl>
<dt>
<code>values</code>
</dt>
<dd>
A named list of RAxML-NG settings.
</dd>
</dl>

## See Also

<a href="https://github.com/amkozlov/raxml-ng/wiki/Input-data">RAxML-NG
input</a>,
<a href="https://github.com/amkozlov/raxml-ng/wiki/Output:-files-and-settings">RAxML-NG
output</a>

## Examples

``` r
library("immundapter")


sequences <- c(
  "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT",
  "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGA",
  "ACGTACGTACGTACGTACGTACGTACGTACGTACGTTCGT",
  "ACGTACGTACGTACGTACGTACGTACGTACGTTCGTACGT"
)
input_file <- tempfile(fileext = ".tsv")
utils::write.table(data.frame(
  cell_id = paste0("cell-", seq_along(sequences)),
  locus = "IGH", sequence = sequences, clone_id = "family-a",
  clone_sequence_alignment = sequences, duplicate_count = 1L
), input_file, sep = "\t", row.names = FALSE, quote = FALSE)
idata <- immundata::read_repertoires(
  path = input_file,
  schema = immundata::make_receptor_schema(
    features = "sequence", chains = "IGH"
  ),
  barcode_col = "cell_id", locus_col = "locus",
  umi_col = "duplicate_count", repertoire_schema = NULL,
  preprocess = NULL, postprocess = NULL, rename_columns = NULL,
  output_folder = tempfile("raxml-immundata-")
)
config <- conf_raxml_ng(
  alignment_col = "clone_sequence_alignment",
  family_col = "clone_id",
  analysis = "search"
)
tool <- tool_raxml_ng(config)
tool$check()

run <- tool$run(idata, output_dir = tempfile("raxml-ng-run-"))
results <- run$results()
dplyr::collect(results$tree_index)
trees <- results$read_trees()
```
