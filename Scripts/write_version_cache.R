# Write a local cache of version and tag history for environments without git
# (e.g. shinyapps.io). Run from the repository root:
#   Rscript Scripts/write_version_cache.R

# Reuse the app's own version logic so the cache always matches what the app
# displays (major.minor aware, e.g. "V1.8").
source("Server/Functions/get_app_version.R")

if (!has_git()) {
  stop("git is not available; run this script locally where git is present.")
}

# Force a fresh tag fetch for the snapshot.
Sys.setenv(APP_GIT_FETCH = "1")
invisible(safe_sys("git", c("fetch", "--tags", "--prune", "--force")))

ordered <- order_version_tags(safe_sys("git", c("tag", "--list")))
version <- if (length(ordered)) version_label(ordered[1]) else "No version tagged"
history <- get_tag_history(limit = 1000)

cache <- list(version = version, history = history)
saveRDS(cache, file = "version_cache.rds")
cat("Wrote version_cache.rds with version:", version, "and", nrow(history), "history rows\n")
