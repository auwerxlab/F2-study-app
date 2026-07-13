has_git <- function() nzchar(Sys.which("git"))

safe_sys <- function(cmd, args = character()) {
  tryCatch(system2(cmd, args, stdout = TRUE, stderr = TRUE), error = function(e) character())
}

# Network fetch of tags is opt-in: deployed apps should not hit the network on
# every startup / Info-tab view. Set APP_GIT_FETCH=1 to refresh tags live.
should_fetch <- function() {
  tolower(Sys.getenv("APP_GIT_FETCH", "false")) %in% c("1", "true", "yes", "on")
}

# Read a locally cached version file if present (used on platforms without git)
read_version_cache <- function(path = "version_cache.rds") {
  tryCatch({
    if (!file.exists(path)) return(NULL)
    readRDS(path)
  }, error = function(e) NULL)
}

# ---- Versioning (major.minor, e.g. "1.8"; legacy integer tags still accepted) ----
# A version tag is "major.minor" ("1.8") or a legacy integer ("7", shown as "1.7").
version_pattern <- "^(v)?([0-9]+\\.[0-9]+|[0-9]+)$"

# Comparable key c(scheme, major, minor): scheme 1 = major.minor (newer than
# any legacy integer), 0 = legacy integer. Used to pick / order versions.
version_key <- function(tag) {
  t <- sub("^v", "", tag, ignore.case = TRUE)
  if (grepl("^[0-9]+\\.[0-9]+$", t)) {
    p <- as.integer(strsplit(t, ".", fixed = TRUE)[[1]])
    return(c(scheme = 1L, major = p[1], minor = p[2]))
  }
  c(scheme = 0L, major = suppressWarnings(as.integer(t)), minor = -1L)
}

# Display label. No "V" prefix - callers add "Version " where they want it. Legacy
# integer tags are shown in the new major.minor scheme: "7" -> "1.7", "1.8" -> "1.8".
# (If you'd rather the first release be 1.0 instead of 1.1, change the map to
#  paste0("1.", as.integer(t) - 1L).)
version_label <- function(tag) {
  t <- sub("^v", "", tag, ignore.case = TRUE)
  if (grepl("^[0-9]+$", t)) t <- paste0("1.", t)   # legacy integer "7" -> "1.7"
  t
}

# Keep only version-like tags, ordered newest-first.
order_version_tags <- function(tags) {
  keep <- grep(version_pattern, tags, value = TRUE, ignore.case = TRUE)
  if (!length(keep)) return(character(0))
  keys <- lapply(keep, version_key)
  ord <- order(vapply(keys, `[[`, integer(1), "scheme"),
               vapply(keys, `[[`, integer(1), "major"),
               vapply(keys, `[[`, integer(1), "minor"),
               decreasing = TRUE)
  keep[ord]
}

# Latest app version (prefers env vars; else local git tags; else cache; e.g. "1.8")
get_app_version <- function() {
  envs <- Sys.getenv(c("APP_VERSION", "CI_COMMIT_TAG", "SOURCE_VERSION"), NA_character_)
  env_version <- envs[!is.na(envs) & nzchar(envs)][1]
  if (!is.na(env_version)) {
    # Normalize "1.8" / "v1.8" / "7" to "1.8" / "1.7"; pass anything else through
    t <- sub("^v", "", env_version, ignore.case = TRUE)
    return(if (grepl("^([0-9]+\\.[0-9]+|[0-9]+)$", t)) version_label(env_version) else env_version)
  }

  # Read local git tags (no network unless explicitly opted in via APP_GIT_FETCH)
  if (has_git()) {
    if (should_fetch()) invisible(safe_sys("git", c("fetch", "--tags", "--prune", "--force")))
    ordered <- order_version_tags(safe_sys("git", c("tag", "--list")))
    if (length(ordered)) return(version_label(ordered[1]))
  }

  # Fall back to the cached version file (e.g. on platforms without git)
  cache <- read_version_cache()
  if (!is.null(cache) && is.character(cache$version) && nzchar(cache$version)) {
    return(cache$version)
  }
  "No version tagged"
}

# Full history: rows = tag (as "1.8"), short sha, date, commit message
get_tag_history <- function(limit = 50) {
  empty <- data.frame(tag = character(), commit = character(), date = character(), message = character())

  # Without git, fall back to the cached history file
  cache_history <- function() {
    cache <- read_version_cache()
    if (!is.null(cache) && is.data.frame(cache$history)) cache$history else empty
  }
  if (!has_git()) return(cache_history())

  # Read local tags (no network unless explicitly opted in via APP_GIT_FETCH)
  if (should_fetch()) invisible(safe_sys("git", c("fetch", "--tags", "--prune", "--force")))
  all_tags <- safe_sys("git", c("tag", "--list"))
  if (!length(all_tags)) return(cache_history())

  # version tags (major.minor or legacy integer), newest first
  keep <- order_version_tags(all_tags)
  if (!length(keep)) return(cache_history())
  if (length(keep) > limit) keep <- keep[seq_len(limit)]

  rows <- lapply(keep, function(tag) {
    ref <- paste0("refs/tags/", tag)
    # resolve to commit (works for lightweight and annotated tags)
    sha <- safe_sys("git", c("rev-list", "-n", "1", ref))
    if (!length(sha)) return(NULL)
    sha <- sha[1]
    msg  <- safe_sys("git", c("log", "-1", "--pretty=%B", sha))
    date <- safe_sys("git", c("log", "-1", "--date=iso", "--pretty=%ad", sha))
    data.frame(
      tag     = version_label(tag),
      commit  = substr(sha, 1, 7),
      date    = if (length(date)) date[1] else "",
      message = if (length(msg)) paste(msg, collapse = "\n") else "",
      stringsAsFactors = FALSE
    )
  })

  rows <- Filter(Negate(is.null), rows)
  if (!length(rows)) return(empty)
  do.call(rbind, rows)
}
