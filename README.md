<div align="center">
  <h1>🔌 <code>immundapter</code> – <strong>Lightweight R Adapters for Immune Repertoire Tools</strong></h1>
</div>

---

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
  <a href="https://immunomind.github.io/immundapter/">Documentation</a>
  |
  <a href="https://immunomind.github.io/docs/">Ecosystem</a>
</p>

---

Analyze `ImmunData` objects or immune receptor files with external tools, read result tables to R, and add
the results as annotations to `ImmunData` objects.

## Supported tools

| Tool | Input | Result | R help |
| --- | --- | --- | --- |
| tcrdist3 | `ImmunData` or a receptor table | T-cell receptor distances and clusters | `?tool_tcrdist3` |
| ImmGLIPH | `ImmunData` with TRB receptors | Groups of similar TRB receptors | `?tool_immgliph` |
| IgBLAST | FASTA with optional annotations | Germline gene calls | `?tool_igblast` |
| RAxML-NG | `ImmunData` with aligned B-cell receptor sequences | Phylogenetic trees by family | `?tool_raxml_ng` |
| VDJdb | `ImmunData` | Antigen specificity annotations | `?data_vdjdb` |

## Install

Use R 4.1 or later.

```r
install.packages("pak")
pak::pkg_install("immunomind/immundapter")
```

## Quick start

1. Print the tcrdist3 setup commands in R.

   ```r
   tool <- immundapter::tool_tcrdist3()
   tool$help()
   ```

2. Run the setup commands in your terminal.

3. Cluster the T-cell receptors in the bundled example.

   ```r
   input <- system.file(
     "extdata", "tcrdist3", "tcr-cluster-simple.tsv",
     package = "immundapter", mustWork = TRUE
   )
   run <- tool$run(input, output_dir = tempdir())
   ```

4. Read the cluster table.

   ```r
   dplyr::collect(run$results()$clusters)
   ```

## Results folder

`immundapter::set_output_dir()` command saves a default folder for tools to read and write results from and to.
