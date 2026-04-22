#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

usage <- function(status = 0L) {
  cat(
    "Usage: Rscript docker/update_dockerhub_overview.R <docker-tag>\n",
    "\n",
    "Example:\n",
    "  Rscript docker/update_dockerhub_overview.R 1.97\n",
    sep = ""
  )
  quit(save = "no", status = status)
}

if (!length(args) || args[1] %in% c("-h", "--help")) {
  usage(0L)
}

if (length(args) != 1L) {
  usage(1L)
}

new_tag <- args[[1]]
if (!grepl("^[0-9]+(\\.[0-9]+)*$", new_tag)) {
  stop("Expected a numeric Docker tag such as '1.96' or '1.96.1'.", call. = FALSE)
}

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
if (!length(script_arg)) {
  stop("Could not determine script location from command arguments.", call. = FALSE)
}

script_path <- normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/", mustWork = TRUE)
docker_dir <- dirname(script_path)
overview_path <- file.path(docker_dir, "DOCKERHUB_OVERVIEW.md")

if (!file.exists(overview_path)) {
  stop("Could not find docker/DOCKERHUB_OVERVIEW.md next to the script.", call. = FALSE)
}

extract_numeric_tags <- function(lines) {
  if (!length(lines)) {
    return(character(0))
  }
  tags <- base::character(0)
  for (line in lines) {
    hit <- regmatches(
      line,
      regexec("^\\s*-\\s*`?([0-9]+(?:\\.[0-9]+)*)`?(?:\\s|$)", line, perl = TRUE)
    )[[1]]
    if (length(hit) >= 2L) {
      tags <- c(tags, hit[[2]])
    }
  }
  tags
}

lines <- readLines(overview_path, warn = FALSE)

recommended_idx <- match("## Recommended tags", lines)
quick_start_idx <- match("## Quick start", lines)

if (is.na(recommended_idx) || is.na(quick_start_idx) || quick_start_idx <= recommended_idx) {
  stop("Could not locate the Docker Hub overview tag section.", call. = FALSE)
}

section_lines <- lines[(recommended_idx + 1L):(quick_start_idx - 1L)]
older_heading_rel <- grep("^Older tags", section_lines)

if (length(older_heading_rel) != 1L) {
  stop("Could not locate the older-tags subsection inside DOCKERHUB_OVERVIEW.md.", call. = FALSE)
}

recommended_lines <- section_lines[seq_len(older_heading_rel - 1L)]
older_lines <- if (older_heading_rel < length(section_lines)) {
  section_lines[(older_heading_rel + 1L):length(section_lines)]
} else {
  character(0)
}

current_recommended <- extract_numeric_tags(recommended_lines)
current_recommended <- if (length(current_recommended)) current_recommended[[1]] else NA_character_

older_tags <- c(
  if (!is.na(current_recommended)) current_recommended else character(0),
  extract_numeric_tags(older_lines)
)
older_tags <- older_tags[!duplicated(older_tags)]
older_tags <- older_tags[older_tags != new_tag]

replacement <- c(
  "## Recommended tags",
  "",
  sprintf("- `%s` for a pinned, reproducible setup", new_tag),
  "- `latest` if you prefer the moving convenience tag",
  "",
  "Older tags still kept for older reproducible runs:",
  ""
)

if (length(older_tags)) {
  replacement <- c(replacement, sprintf("- `%s`", older_tags))
}

replacement <- c(replacement, "")

updated_lines <- c(
  if (recommended_idx > 1L) lines[seq_len(recommended_idx - 1L)] else character(0),
  replacement,
  lines[quick_start_idx:length(lines)]
)

updated_lines <- gsub(
  "docker pull therealtomek/hcocena:[0-9]+(\\.[0-9]+)*",
  sprintf("docker pull therealtomek/hcocena:%s", new_tag),
  updated_lines
)
updated_lines <- gsub(
  "docker run --rm -p 8787:8787 -e PASSWORD=hcocena therealtomek/hcocena:[0-9]+(\\.[0-9]+)*",
  sprintf("docker run --rm -p 8787:8787 -e PASSWORD=hcocena therealtomek/hcocena:%s", new_tag),
  updated_lines
)

writeLines(updated_lines, overview_path, useBytes = TRUE)

cat("Updated:", overview_path, "\n", sep = "")
cat("Recommended Docker tag:", new_tag, "\n", sep = " ")
if (!is.na(current_recommended) && current_recommended != new_tag) {
  cat("Previous pinned tag moved to older tags:", current_recommended, "\n", sep = " ")
} else if (!is.na(current_recommended) && current_recommended == new_tag) {
  cat("Pinned tag already matched requested tag before the update.\n")
}
if (length(older_tags)) {
  cat("Older reproducibility tags:", paste(older_tags, collapse = ", "), "\n", sep = " ")
}
cat(
  "Next step: after pushing the new Docker image tag, sync the contents of ",
  "docker/DOCKERHUB_OVERVIEW.md to the Docker Hub Overview page if it is not auto-synced.\n",
  sep = ""
)
