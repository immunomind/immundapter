make_test_tool <- function(name = "test-tool") {
  TestToolConfig <- R6::R6Class(
    classname = paste0("TestToolConfig", sample.int(1000000L, 1L)),
    inherit = AbstractConfig,
    cloneable = FALSE,
    public = list(
      value = NULL,
      initialize = function(value = "default") {
        checkmate::assert_string(value)
        super$initialize(list(value = value))
      }
    )
  )

  config_class <- TestToolConfig
  TestToolRun <- R6::R6Class(
    classname = paste0("TestToolRun", sample.int(1000000L, 1L)),
    inherit = ToolRun,
    cloneable = FALSE,
    private = list(
      results_impl = function() {
        ToolResults$new(self)
      }
    )
  )
  run_class <- TestToolRun
  TestTool <- R6::R6Class(
    classname = paste0("TestTool", sample.int(1000000L, 1L)),
    inherit = AbstractTool,
    cloneable = FALSE,
    public = list(
      initialize = function(config = config_class$new()) {
        super$initialize(name, config, run_class)
      }
    ),
    private = list(
      help_impl = function() {
        "test help"
      },
      validate_input_impl = function(input, config) {
        if (inherits(input, "ImmunData")) {
          return(input)
        }
        checkmate::assert_string(input, min.chars = 1L)
        checkmate::assert_file_exists(input)
        fs::path_abs(input)
      },
      run_impl = function(input, config, run_dir) {
        list(run_dir = run_dir)
      }
    )
  )

  list(
    tool = TestTool$new(),
    config_class = TestToolConfig
  )
}

test_that("AbstractTool exposes help instead of install", {
  tool <- make_test_tool()$tool

  expect_identical(tool$help(), "test help")
  expect_true("help" %in% names(tool))
  expect_false("install" %in% names(tool))
})

test_that("AbstractTool runs with a file and explicit output directory", {
  root <- withr::local_tempdir()
  withr::local_envvar(R_USER_CONFIG_DIR = root)
  input <- fs::path(root, "input.tsv")
  writeLines("value", input)
  set_output_dir(fs::path(root, "default"))
  output_dir <- fs::path(root, "results")
  fixture <- make_test_tool()

  run <- fixture$tool$run(input, output_dir = output_dir)

  expect_true(inherits(run, "ToolRun"))
  expect_equal(run$status, "success")
  expect_equal(run$input$type, "files")
  expect_equal(unname(run$input$paths), as.character(fs::path_abs(input)))
  expect_equal(
    as.character(fs::path_dir(run$path)),
    as.character(fs::path_abs(output_dir))
  )
  expect_true(fs::dir_exists(run$path))
  expect_error(
    run$status <- "partial",
    "`status` is read-only.",
    fixed = TRUE
  )
  expect_true(inherits(run$results(), "ToolResults"))
})

test_that("ToolResults always exposes an informative annotate method", {
  root <- withr::local_tempdir()
  run <- ToolRun$new(list(
    status = "success",
    tool = "test-tool",
    path = root,
    artifacts = list(),
    config = list(),
    input = list(type = "object")
  ))
  TestImmunData <- R6::R6Class(
    classname = paste0("ImmunData", sample.int(1000000L, 1L))
  )
  idata <- TestImmunData$new()
  class(idata) <- unique(c("ImmunData", class(idata)))
  results <- ToolResults$new(run)

  expect_error(
    run$results(),
    "must implement `$results()` for this run type",
    fixed = TRUE
  )
  expect_true("annotate" %in% names(results))
  expect_error(
    results$annotate(idata),
    "cannot annotate `ImmunData`",
    class = "immundapter_unannotatable_results"
  )
})

test_that("AbstractTool requires a destination for file runs without a preference", {
  root <- withr::local_tempdir()
  withr::local_envvar(R_USER_CONFIG_DIR = root)
  input <- fs::path(root, "input.tsv")
  writeLines("value", input)
  fixture <- make_test_tool()

  expect_error(fixture$tool$run(input), "set_output_dir()", fixed = TRUE)
  expect_error(fixture$tool$run(input, output_dir = ""), "output_dir")
  expect_false(fs::dir_exists(tools::R_user_dir("immundapter", "config")))
})

test_that("AbstractTool uses the saved absolute directory for file runs", {
  root <- withr::local_tempdir()
  withr::local_envvar(c(R_USER_CONFIG_DIR = root, R_FS_HOME = root))
  input <- fs::path(root, "input.tsv")
  writeLines("value", input)
  fixture <- make_test_tool()
  destination <- fs::path(root, "default runs")

  withr::with_dir(root, set_output_dir("default runs"))
  expect_false(fs::dir_exists(destination))

  run <- fixture$tool$run(input)

  expect_equal(fs::path_real(fs::path_dir(run$path)), fs::path_real(destination))
  expect_true(fs::dir_exists(run$path))

  new_destination <- fs::path(root, "new runs")
  expect_equal(set_output_dir("~/new runs"), new_destination)
  run <- fixture$tool$run(input)
  expect_equal(fs::path_real(fs::path_dir(run$path)), fs::path_real(new_destination))
})

test_that("AbstractTool uses the ImmunData snapshot artifacts path", {
  root <- withr::local_tempdir()
  withr::local_envvar(R_USER_CONFIG_DIR = root)
  set_output_dir(fs::path(root, "default"))
  artifacts_path <- fs::path(root, "artifacts", "analysis", "v001")

  TestImmunData <- R6::R6Class(
    classname = "ImmunData",
    public = list(
      provenance = NULL,
      initialize = function(provenance) {
        self$provenance <- provenance
      }
    )
  )
  input <- TestImmunData$new(list(
    current_path = fs::path(root, "snapshots", "analysis", "v001"),
    artifacts_path = artifacts_path,
    snapshot_id = "snapshot-1"
  ))
  fixture <- make_test_tool()

  run <- fixture$tool$run(input)

  expect_equal(
    as.character(fs::path_dir(run$path)),
    as.character(fs::path_abs(artifacts_path))
  )
  expect_equal(run$input$type, "immundata")
  expect_equal(run$input$snapshot_id, "snapshot-1")
})

test_that("AbstractTool requires output_dir for unpersisted ImmunData", {
  root <- withr::local_tempdir()
  withr::local_envvar(R_USER_CONFIG_DIR = root)
  set_output_dir(fs::path(root, "default"))
  TestImmunData <- R6::R6Class(
    classname = "ImmunData",
    public = list(
      provenance = list(
        current_path = NULL,
        artifacts_path = NULL,
        snapshot_id = NULL
      )
    )
  )
  fixture <- make_test_tool()

  expect_error(
    fixture$tool$run(TestImmunData$new()),
    "has not been saved to disk",
    fixed = TRUE
  )
})

test_that("AbstractTool uses an immutable per-run configuration", {
  root <- withr::local_tempdir()
  input <- fs::path(root, "input.tsv")
  writeLines("value", input)
  fixture <- make_test_tool()
  config <- fixture$config_class$new("custom")

  run <- fixture$tool$run(
    input,
    config = config,
    output_dir = fs::path(root, "results")
  )

  expect_equal(run$config$value, "custom")
  expect_error(
    config$value <- "changed-after-run",
    "locked binding",
    fixed = TRUE
  )
})
