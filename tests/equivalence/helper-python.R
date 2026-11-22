# Python references, run in a separate process. Nothing here calls hypergraphs.
#
# The interpreter is taken from the environment: HYPERGRAPHS_REFERENCE_PYTHON (the verify job sets it to the
# virtualenv that holds the reference modules on the ci machine), else the helper's verify_python()
# (VERIFY_PYTHON, RETICULATE_PYTHON, ~/equiv-py, then `python3` on the PATH), so a re-run in Prova finds the
# same virtualenv on a machine set up like ci. A module that the interpreter cannot import is a recorded skip,
# named `python:<module>`.

#' The Python interpreter of the references
#'
#' @return A length-one character vector, empty when no interpreter is found.
ref_python <- function() {
  given <- Sys.getenv("HYPERGRAPHS_REFERENCE_PYTHON", unset = "")
  if (nzchar(given)) return(path.expand(given))
  verify_python()
}

#' Whether the interpreter imports a module
#'
#' @param module The import name.
#' @return `TRUE` or `FALSE`.
ref_python_has <- function(module) {
  python <- ref_python()
  nzchar(python) && file.exists(python) &&
    identical(as.integer(system2(python, c("-c", shQuote(sprintf("import %s", module))),
                                 stdout = FALSE, stderr = FALSE)), 0L)
}

#' Skip a test whose Python reference cannot be imported, and record the skip
#'
#' The helper of the standard records skips for R packages only (`skip_without_reference()` asks
#' `requireNamespace()`); a module that is absent is handed to it as `python:<name>`, which no R package is
#' called, so the skip is recorded under that name.
#'
#' @param modules A named character vector: names are the references (as in `verify.json`), values the import
#'   names (`c("scikit-learn" = "sklearn")`); an unnamed value is both.
#' @return `TRUE` invisibly when every module imports.
skip_without_python <- function(modules) {
  stopifnot("`modules` must be module names" = is.character(modules) && length(modules) >= 1L)
  references <- names(modules) %||% modules
  references[!nzchar(references)] <- modules[!nzchar(references)]
  absent <- references[!vapply(modules, ref_python_has, logical(1))]
  if (length(absent)) skip_without_reference(paste0("python:", absent))
  invisible(TRUE)
}

#' The installed version of a Python module, for the note of a comparison
#'
#' @param distribution The distribution name (`"scikit-learn"`, `"xgi"`).
#' @return A length-one character vector, e.g. `"xgi 1.1"`.
ref_python_version <- function(distribution) {
  out <- system2(ref_python(), c("-c", shQuote(sprintf("import importlib.metadata as m; print(m.version('%s'))", distribution))),
                 stdout = TRUE, stderr = TRUE)
  sprintf("Python %s %s", distribution, utils::tail(out, 1L))
}

#' Run a Python script on an input and read back what it writes
#'
#' The input is written as JSON to a file whose path is the script's first argument; the script writes its
#' result as JSON to the path in its second argument. A script that fails stops the test with its output.
#'
#' @param script The Python code, a character vector of lines.
#' @param input A list, written with `jsonlite::write_json(auto_unbox = TRUE, digits = NA)`.
#' @return The parsed result (`jsonlite::fromJSON(simplifyVector = TRUE)`).
ref_python_run <- function(script, input) {
  files <- c(script = tempfile("reference_", fileext = ".py"), input = tempfile("input_", fileext = ".json"),
             output = tempfile("output_", fileext = ".json"))
  on.exit(unlink(files), add = TRUE, after = FALSE)
  writeLines(script, files[["script"]])
  jsonlite::write_json(input, files[["input"]], auto_unbox = TRUE, digits = NA)
  out <- system2(ref_python(), c("-I", shQuote(files[["script"]]), shQuote(files[["input"]]), shQuote(files[["output"]])),
                 stdout = TRUE, stderr = TRUE)
  status <- attr(out, "status") %||% 0L
  if (!identical(as.integer(status), 0L) || !file.exists(files[["output"]])) {
    stop(errorCondition(paste(c("The Python reference failed:", out), collapse = "\n"),
                        class = "reference_python_failed", call = NULL))
  }
  # Python's json writes NaN and Infinity, which are not JSON: NaN is read as a missing value, an infinity as one
  text <- paste(readLines(files[["output"]], warn = FALSE), collapse = "\n")
  text <- gsub("-Infinity\\b", "-1e999", gsub("(?<![-\\w])Infinity\\b", "1e999", text, perl = TRUE))
  jsonlite::fromJSON(gsub("\\bNaN\\b", "null", text), simplifyVector = TRUE)
}
