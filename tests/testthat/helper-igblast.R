igblast_test_submitted_queries <- function() {
  data.frame(
    igblast_query_id = c("query-1", "query-2", "query-3"),
    sequence = c("ACGT", "TGCA", "AAAA")
  )
}


read_igblast_fixture <- function() {
  path <- system.file(
    "extdata",
    "igblast",
    "human-trb.tsv",
    package = "immundapter"
  )
  stopifnot(nzchar(path))

  duckplyr::read_csv_duckdb(
    path,
    options = list(delim = "\t", header = TRUE)
  )
}


igblast_test_runtime_info <- function(
  germline_db = "IMGT.mock.Homo_sapiens.TRB",
  c_region_db = ""
) {
  list(
    package_installed = TRUE,
    igblastr_version = "1.0.23",
    igblast_installed = TRUE,
    igblast_version = "1.22.0",
    igblast_build = "igblast 1.22.0, mock build",
    igblast_root = "/mock/igblast",
    germline_db = germline_db,
    c_region_db = c_region_db
  )
}


read_igblast_test_queries <- function(arguments) {
  fasta <- readLines(arguments$query)

  data.frame(
    sequence_id = sub(
      "^>",
      "",
      fasta[seq.int(1L, length(fasta), by = 2L)]
    ),
    sequence = fasta[seq.int(2L, length(fasta), by = 2L)]
  )
}


write_igblast_test_airr <- function(arguments, airr) {
  utils::write.table(
    airr,
    file = arguments$out,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = ""
  )
  invisible(NULL)
}
