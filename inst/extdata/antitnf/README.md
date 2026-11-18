# GSE290006 anti-TNF vaccine-response BCR dataset

This is a compact, biologically interpretable single-cell BCR dataset for the
immundapter end-to-end tutorials. It contains spike-specific memory B cells
collected 3–4 months after SARS-CoV-2 mRNA vaccine doses 2 and 3 from:

- 3 healthy controls;
- 4 people with inflammatory bowel disease (IBD) receiving anti-IL-12/23;
- 5 people with IBD receiving anti-TNF.

The same participants were sampled at both timepoints except `TNF-05`, who was
sampled only after dose 2. Keeping all three groups makes it possible to ask the
study's real biological question: how anti-TNF therapy changes antigen-specific
germinal-center output, using anti-IL-12/23-treated IBD and healthy participants
as meaningful comparators.

## Files

- `bcr_contigs.fasta.gz`: nucleotide contigs for rerunning IgBLAST.
- `contig_annotations.csv.gz`: combined Cell Ranger V(D)J annotations with
  participant, treatment group, and timepoint.
- `igblast_airr.tsv.gz`: the study's IgBLAST/Change-O AIRR table, with the same
  participant metadata, for workflows that start after IgBLAST.
- `cells.tsv.gz`: one row per retained B cell with HTO assignment, GEX/HTO QC,
  productive chain counts, and an exact heavy-light pairing flag.
- `participants.tsv`: hashtag-to-participant and treatment-group mapping.

All cell and contig identifiers are prefixed with `T4_` or `T6_` to make them
globally unique after combining the two 10x runs. `T4` is 3–4 months after dose
2; `T6` is 3–4 months after dose 3. Cell Ranger clonotype identifiers are also
timepoint-prefixed because those calls are local to a 10x run.

The FASTA and Cell Ranger table include every contig from a cell retained by
the BCR and hashtag filters, including nonproductive contigs. The AIRR file is
the corresponding IgBLAST `db-pass` output and can therefore contain fewer
records. `exact_heavy_light_pair` means exactly one productive, high-confidence
IGH and exactly one productive, high-confidence IGK or IGL contig.

## Selection and provenance

Source: [NCBI GEO GSE290006](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE290006),
specifically BCR samples `GSM8803450` (T4) and `GSM8803453` (T6), plus the T4
and T6 GEX/HTO filtered feature matrices used to recover participant identity.

The HTO demultiplexing reproduces the authors' public analysis: RNA QC of
`200 < nFeature_RNA < 6500` and mitochondrial reads below 10%, HTO CLR
normalization across cells, `HTODemux(positive.quantile = 0.999, seed = 25)`,
and `HTO_margin >= 2`. Cells are retained here when they also have at least one
productive, high-confidence BCR contig. The complete preparation is in
`data-raw/prepare-gse290006.R`.

This produces 1,364 BCR cells, of which 1,223 have an exact productive
heavy-light pair. These are real antigen-specific repertoires rather than a
synthetic phylogeny fixture. Individual expanded lineages are consequently
small; lineage inference should be performed within participant, and selected
expanded or cross-timepoint lineages are appropriate for the RAxML-NG section.

## Attribution

Cheung MW et al. *Reduced vaccine-induced germinal center outputs in patients
with inflammatory bowel disease treated with anti-TNF biologics.* Journal of
Clinical Investigation (2025). [doi:10.1172/JCI192589](https://doi.org/10.1172/JCI192589).

The authors' analysis code is available in
[Cheung_JCI2025](https://github.com/michelle-cheung25/Cheung_JCI2025).
