args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2L) {
  stop("usage: Rscript install_blood_composition_dependencies.R <project-library> <audit-json>")
}

project_library <- normalizePath(args[[1]], winslash = "/", mustWork = FALSE)
audit_path <- args[[2]]
dir.create(project_library, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(project_library, .libPaths()))

methods <- list(
  xCell = list(
    package = "xCell",
    version = "1.1.0",
    repository = "https://github.com/dviraran/xCell",
    commit = "20e2919eefd37e15af35f29f4944e30697098a28",
    install_subdir = "."
  ),
  MCPcounter = list(
    package = "MCPcounter",
    version = "1.2.0",
    repository = "https://github.com/ebecht/MCPcounter",
    commit = "b6eac73e91c246fcff0bb1a5c68a816cd588fc48",
    install_subdir = "Source"
  )
)

json_escape <- function(value) {
  value <- gsub("\\\\", "\\\\\\\\", value)
  value <- gsub("\"", "\\\\\"", value)
  value <- gsub("\r", "\\\\r", value)
  gsub("\n", "\\\\n", value)
}

json_string <- function(value) {
  paste0("\"", json_escape(as.character(value)), "\"")
}

json_bool <- function(value) {
  if (isTRUE(value)) "true" else "false"
}

json_object <- function(fields) {
  paste0("{", paste(paste0(json_string(names(fields)), ":", unlist(fields)), collapse = ","), "}")
}

write_audit <- function(method_results) {
  rendered_methods <- vapply(
    names(method_results),
    function(name) json_object(method_results[[name]]),
    character(1)
  )
  audit <- json_object(list(
    status = json_string("INSTALL_AUDIT_ONLY"),
    timestamp_utc = json_string(format(Sys.time(), tz = "UTC", usetz = TRUE)),
    r_version = json_string(R.version.string),
    r_library_scope = json_string("project-local-only"),
    project_library = json_string(project_library),
    attempted_without_system_library_write = json_bool(TRUE),
    methods = paste0("{", paste(paste0(json_string(names(rendered_methods)), ":", rendered_methods), collapse = ","), "}"),
    session_info = json_string(paste(capture.output(sessionInfo()), collapse = "\n"))
  ))
  writeLines(audit, audit_path, useBytes = TRUE)
}

method_results <- lapply(methods, function(method) {
  list(
    package = json_string(method$package),
    requested_version = json_string(method$version),
    repository = json_string(method$repository),
    resolved_commit = json_string(method$commit),
    install_subdir = json_string(method$install_subdir),
    status = json_string("BLOCKED"),
    reason = json_string("not attempted")
  )
})

if (!requireNamespace("remotes", quietly = TRUE)) {
  for (name in names(method_results)) {
    method_results[[name]]$reason <- json_string("project-local remotes package unavailable; no system-library write attempted")
  }
  write_audit(method_results)
  quit(status = 0L)
}

for (name in names(methods)) {
  method <- methods[[name]]
  result <- tryCatch(
    {
      remotes::install_github(
        repo = sub("https://github.com/", "", method$repository),
        ref = method$commit,
        subdir = method$install_subdir,
        lib = project_library,
        upgrade = "never",
        dependencies = TRUE,
        quiet = TRUE
      )
      installed <- requireNamespace(method$package, quietly = TRUE, lib.loc = project_library)
      if (!installed) stop("installation returned without the requested package")
      list(status = "READY", reason = "fixed commit installed in project-local R library")
    },
    error = function(error) list(status = "BLOCKED", reason = conditionMessage(error))
  )
  method_results[[name]]$status <- json_string(result$status)
  method_results[[name]]$reason <- json_string(result$reason)
}

write_audit(method_results)
