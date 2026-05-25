#' Time coarsening and follow-up augmentation for start-stop survival data
#'
#' @description
#' `tcoarsen()` harmonizes irregular longitudinal survival data in start-stop
#' form onto a user-specified discrete time grid. It maps observed update times
#' to grid landmarks using either floor or ceiling coarsening, rebuilds start-stop
#' intervals, and augments each individual's follow-up record by carrying the
#' last known covariate and treatment values forward to empty grid intervals.
#'
#' The function is intended as a transparent preprocessing utility for real-world
#' clinical trial datasets with sparse or irregular follow-up times.The grid width
#' and coarsening direction are substantive analysis choices and should usually be
#' prespecified and examined in sensitivity analysis.
#'
#' For use before [multimsm()], include the triplet variables such as `rand`,
#' `cross`, and `subseq` in `covs`, and set `absorb_vars = c("cross", "subseq")`
#' when crossover and subsequent therapy indicators should be interpreted as
#' absorbing current-status variables. The output remains an ordinary data frame
#' that can be passed directly to downstream modeling functions.
#'
#' @details
#' Suppose an individual has observed rows indexed by `j`, representing intervals
#' `[start_ij, stop_ij)` over which the row values are assumed to apply. Given
#' grid origin `origin` and bin width `bin_width`, an observed update time `s` is
#' mapped to
#'
#' * `origin + floor((s - origin) / bin_width) * bin_width` when
#'   `dir_coarsen = "floor"`;
#' * `origin + ceiling((s - origin) / bin_width) * bin_width` when
#'   `dir_coarsen = "ceiling"`.
#'
#' After snapping update times, `tcoarsen()` inserts one row for each missing grid
#' interval during which the subject remains under observation. Values of all
#' carried variables are filled by last-known-value-carried-forward (LKCF).
#' Terminal event/censoring times are preserved exactly by default, so an event
#' at day 173 can remain `stop = 173` even if the coarsening grid is monthly.
#'
#' If multiple observed rows map to the same coarsened update time, the last
#' observed row within that coarsened bin is kept. This rule is simple and
#' reproducible, but it necessarily discards within-bin ordering information.
#'
#' @param data A `data.frame` in start-stop long format, with one or more rows
#'   per subject.
#' @param id A character string for the subject identifier variable name.
#' @param start A character string naming the interval start-time
#'   variable.
#' @param stop A character string for the interval stop-time variable name.
#' @param event A character string for the terminal event indicator.
#'   The event variable must be coded `0`/`1`, with at most one event per
#'   subject.
#' @param covs Optional character vector of treatment or covariate variables
#'   to document as carried-forward predictors. Both baseline-only and
#'   time-varying variables are allowed. The function retains all original data
#'   columns, but `covs` is checked for existence and stored in the returned
#'   object's attributes for reproducibility. In practice, all __time-varying__
#'   covariates are encouraged to be included here.
#' @param bin_width A positive numeric value giving the grid width in the
#'   same time units as `start` and `stop`, for example `30` for a 30-day grid.
#' @param dir_coarsen Coarsening direction, either `"floor"` or `"ceiling"`.
#'   Floor coarsening moves updates backward to the start of the containing bin;
#'   ceiling coarsening moves updates forward to the end of the containing bin.
#' @param origin Optional numeric grid origin. If `NULL`, the minimum observed
#'   start time is used. In trial applications this is often `0`.
#' @param absorb_vars Optional character vector of variables to convert to
#'   absorbing status using within-subject cumulative maxima after sorting by
#'   time. These variables must be coded `0`/`1`. For multi-way switching data,
#'   this will often be `c("cross", "subseq")`.
#' @param keep_terminal_time Logical. If `TRUE` (default), each subject's exact
#'   event/censoring exit time is preserved as the final stop time. If `FALSE`,
#'   the exit time is also snapped to the selected grid. The default is strongly
#'   recommended for survival analyses.
#' @param gap_action How to handle gaps between adjacent within-subject
#'   start-stop intervals before coarsening. One of `"stop"`, `"warn"`, or
#'   `"ignore"`. Overlapping intervals always stop with an error.
#' @param add_visit Logical. If `TRUE` (default), add a discrete grid index to
#'   the output.
#' @param visit_name A single character string giving the name of the visit-index
#'   variable to add when `add_visit = TRUE`. Defaults to `"visit"`.
#' @param diagnostics Logical. If `TRUE` (default), attach preprocessing metadata
#'   and a by-visit row/event summary in `diagnostics` object of output.
#' @param verbose Logical. If `TRUE`, print a short preprocessing message.
#'
#' @return
#' A `tcoarsen` object with two components. The first component is a `data.frame`
#' called `dat_coarsen` the time-coarsened input data with updated interval times
#' and time-varying covariates under coarsened time grid; the second component is
#' `diagnostics` for diagnostics settings and basic diagnostics.
#'
#' @seealso [multimsm()]
#'
#' @references
#' Guerra SF, Schnitzer ME, Forget A, Blais L. Impact of discretization of the
#' timeline for longitudinal causal inference methods. \emph{Statistics in
#' Medicine}. 2020 Sept;39(27):0277-6715.
#'
#' @examples
#' if (requireNamespace("survival", quietly = TRUE)) {
#'   data("heart", package = "survival")
#'
#'   heart_30 <- tcoarsen(
#'     data = heart,
#'     id = "id",
#'     start = "start",
#'     stop = "stop",
#'     event = "event",
#'     covs = c("transplant", "age", "surgery"),
#'     bin_width = 30,
#'     dir_coarsen = "floor",
#'     origin = 0,
#'     absorb_vars = "transplant",
#'     verbose = FALSE
#'   )
#'
#'   utils::head(heart_30$dat_coarsen, 10)
#' }
#'
#' @export

tcoarsen <- function(data,
                     id,
                     start,
                     stop,
                     event,
                     covs = NULL,
                     bin_width,
                     dir_coarsen = c("floor", "ceiling"),
                     origin = NULL,
                     absorb_vars = NULL,
                     keep_terminal_time = TRUE,
                     gap_action = c("stop", "warn", "ignore"),
                     add_visit = TRUE,
                     visit_name = "visit",
                     diagnostics = TRUE,
                     verbose = TRUE) {

  dir_coarsen <- match.arg(dir_coarsen)
  gap_action <- match.arg(gap_action)

  if (!is.data.frame(data)) {
    stop("`data` must be a data.frame.", call. = FALSE)
  }

  .tcoarsen_check_string(id, "id")
  .tcoarsen_check_string(start, "start")
  .tcoarsen_check_string(stop, "stop")
  .tcoarsen_check_string(event, "event")
  .tcoarsen_check_string(visit_name, "visit_name")

  needed <- c(id, start, stop, event)
  missing_needed <- setdiff(needed, names(data))
  if (length(missing_needed) > 0L) {
    stop("The following required variables are not in `data`: ",
         paste(missing_needed, collapse = ", "), call. = FALSE)
  }

  if (!is.null(covs)) {
    if (!is.character(covs)) {
      stop("`covs` must be NULL or a character vector of variable names.",
           call. = FALSE)
    }
    missing_covs <- setdiff(covs, names(data))
    if (length(missing_covs) > 0L) {
      stop("The following variables in `covs` are not in `data`: ",
           paste(missing_covs, collapse = ", "), call. = FALSE)
    }
  }

  if (!is.null(absorb_vars)) {
    if (!is.character(absorb_vars)) {
      stop("`absorb_vars` must be NULL or a character vector of variable names.",
           call. = FALSE)
    }
    missing_abs <- setdiff(absorb_vars, names(data))
    if (length(missing_abs) > 0L) {
      stop("The following variables in `absorb_vars` are not in `data`: ",
           paste(missing_abs, collapse = ", "), call. = FALSE)
    }
    absorb_vars <- unique(absorb_vars)
  }

  if (!is.numeric(bin_width) || length(bin_width) != 1L ||
      !is.finite(bin_width) || bin_width <= 0) {
    stop("`bin_width` must be a single positive numeric value.", call. = FALSE)
  }

  if (!is.logical(keep_terminal_time) || length(keep_terminal_time) != 1L ||
      is.na(keep_terminal_time)) {
    stop("`keep_terminal_time` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(add_visit) || length(add_visit) != 1L || is.na(add_visit)) {
    stop("`add_visit` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(diagnostics) || length(diagnostics) != 1L || is.na(diagnostics)) {
    stop("`diagnostics` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("`verbose` must be TRUE or FALSE.", call. = FALSE)
  }

  if (isTRUE(add_visit) && visit_name %in% names(data) &&
      !identical(visit_name, start) && !identical(visit_name, stop) &&
      !identical(visit_name, event)) {
    warning("`visit_name = \"", visit_name, "\"` already exists in `data`; ",
            "it will be overwritten in the returned data.", call. = FALSE)
  }

  df <- data
  n_original <- nrow(df)
  if (n_original == 0L) {
    stop("`data` has zero rows.", call. = FALSE)
  }

  df$.tcoarsen_orig_row <- seq_len(nrow(df))
  df[[start]] <- .tcoarsen_to_numeric_time(df[[start]], start)
  df[[stop]] <- .tcoarsen_to_numeric_time(df[[stop]], stop)
  df[[event]] <- .tcoarsen_to_binary01(df[[event]], event, allow_na = FALSE)

  if (anyNA(df[[id]])) {
    stop("`id` cannot contain missing values.", call. = FALSE)
  }
  if (any(!is.finite(df[[start]]) | !is.finite(df[[stop]]))) {
    stop("`start` and `stop` must contain finite numeric values.", call. = FALSE)
  }
  if (any(df[[stop]] <= df[[start]])) {
    stop("Found rows with `stop <= start`.", call. = FALSE)
  }

  if (is.null(origin)) {
    origin <- min(df[[start]], na.rm = TRUE)
  }
  if (!is.numeric(origin) || length(origin) != 1L || !is.finite(origin)) {
    stop("`origin` must be NULL or a single finite numeric value.", call. = FALSE)
  }

  if (any(df[[start]] < origin)) {
    stop("Found interval start times smaller than `origin`. Choose an earlier ",
         "origin or leave `origin = NULL`.", call. = FALSE)
  }

  for (v in absorb_vars) {
    df[[v]] <- .tcoarsen_to_binary01(df[[v]], v, allow_na = FALSE)
  }

  df <- df[order(df[[id]], df[[start]], df[[stop]], df$.tcoarsen_orig_row),
           , drop = FALSE]

  ev_by_id <- tapply(df[[event]], df[[id]], function(x) sum(x, na.rm = TRUE))
  if (any(ev_by_id > 1L, na.rm = TRUE)) {
    stop("More than one terminal event detected for at least one subject.",
         call. = FALSE)
  }

  .tcoarsen_check_intervals(df = df,
                            id = id,
                            start = start,
                            stop = stop,
                            gap_action = gap_action,
                            tol = sqrt(.Machine$double.eps))

  if (length(absorb_vars) > 0L) {
    for (v in absorb_vars) {
      df[[v]] <- stats::ave(as.integer(df[[v]]), df[[id]], FUN = function(x) cummax(x))
    }
  }

  snap_time <- switch(
    dir_coarsen,
    floor = function(t) origin + floor((t - origin) / bin_width) * bin_width,
    ceiling = function(t) origin + ceiling((t - origin) / bin_width) * bin_width
  )

  diag_env <- new.env(parent = emptyenv())
  diag_env$n_collapsed_rows <- 0L
  diag_env$n_collapsed_bins <- 0L
  diag_env$n_dropped_updates_after_exit <- 0L
  diag_env$n_subjects_with_post_exit_updates <- 0L
  diag_env$n_records_after_bin_collapse <- 0L

  split_ids <- split(df, df[[id]], drop = TRUE)

  out_list <- lapply(split_ids, function(d) {
    d <- d[order(d[[start]], d[[stop]], d$.tcoarsen_orig_row), , drop = FALSE]

    entry_time <- min(d[[start]], na.rm = TRUE)
    final_observed_stop <- max(d[[stop]], na.rm = TRUE)
    has_event <- any(d[[event]] == 1L, na.rm = TRUE)
    event_time <- if (has_event) {
      min(d[[stop]][d[[event]] == 1L], na.rm = TRUE)
    } else {
      Inf
    }

    exit_raw <- if (has_event) event_time else final_observed_stop
    exit_time <- if (isTRUE(keep_terminal_time)) {
      exit_raw
    } else {
      z <- snap_time(exit_raw)
      if (z <= entry_time) {
        z <- entry_time + bin_width
      }
      z
    }

    d$.tcoarsen_snap_start <- snap_time(d[[start]])
    d$.tcoarsen_snap_start[1L] <- entry_time
    d$.tcoarsen_snap_start <- pmax(d$.tcoarsen_snap_start, entry_time)

    drop_after_exit <- d$.tcoarsen_snap_start >= exit_time
    if (any(drop_after_exit)) {
      diag_env$n_dropped_updates_after_exit <-
        diag_env$n_dropped_updates_after_exit + sum(drop_after_exit)
      diag_env$n_subjects_with_post_exit_updates <-
        diag_env$n_subjects_with_post_exit_updates + 1L
    }
    d <- d[!drop_after_exit, , drop = FALSE]
    if (nrow(d) == 0L) {
      return(NULL)
    }

    d <- d[order(d$.tcoarsen_snap_start,
                 d[[start]],
                 d[[stop]],
                 d$.tcoarsen_orig_row), , drop = FALSE]

    duplicated_start <- duplicated(d$.tcoarsen_snap_start) |
      duplicated(d$.tcoarsen_snap_start, fromLast = TRUE)
    if (any(duplicated_start)) {
      collapsed_bins <- length(unique(d$.tcoarsen_snap_start[duplicated_start]))
      n_before <- nrow(d)
      keep_last <- !duplicated(d$.tcoarsen_snap_start, fromLast = TRUE)
      d <- d[keep_last, , drop = FALSE]
      diag_env$n_collapsed_bins <- diag_env$n_collapsed_bins + collapsed_bins
      diag_env$n_collapsed_rows <- diag_env$n_collapsed_rows + (n_before - nrow(d))
    }

    d <- d[order(d$.tcoarsen_snap_start, d$.tcoarsen_orig_row), , drop = FALSE]

    if (length(absorb_vars) > 0L) {
      for (v in absorb_vars) {
        d[[v]] <- cummax(as.integer(d[[v]]))
      }
    }

    diag_env$n_records_after_bin_collapse <-
      diag_env$n_records_after_bin_collapse + nrow(d)

    k_min <- ceiling((entry_time - origin) / bin_width)
    k_max <- floor((exit_time - origin) / bin_width)
    if (k_max >= k_min) {
      grid_starts <- origin + seq.int(from = k_min, to = k_max) * bin_width
      grid_starts <- grid_starts[grid_starts >= entry_time & grid_starts < exit_time]
    } else {
      grid_starts <- numeric(0L)
    }

    augmented_starts <- sort(unique(c(entry_time,
                                      grid_starts,
                                      d$.tcoarsen_snap_start)))
    augmented_starts <- augmented_starts[augmented_starts < exit_time]
    if (length(augmented_starts) == 0L) {
      return(NULL)
    }

    record_times <- d$.tcoarsen_snap_start
    idx <- findInterval(augmented_starts, record_times)
    idx[idx < 1L] <- 1L

    out <- d[idx, , drop = FALSE]
    out[[start]] <- augmented_starts
    out[[stop]] <- c(augmented_starts[-1L], exit_time)
    out[[event]] <- 0L

    if (has_event) {
      out[[event]][nrow(out)] <- 1L
      if (isTRUE(keep_terminal_time)) {
        out[[stop]][nrow(out)] <- event_time
      }
    }

    out <- out[out[[stop]] > out[[start]], , drop = FALSE]
    if (nrow(out) == 0L) {
      return(NULL)
    }

    out
  })

  out_list <- Filter(Negate(is.null), out_list)
  if (length(out_list) == 0L) {
    stop("No valid follow-up intervals remain after coarsening.", call. = FALSE)
  }

  out <- do.call(rbind, out_list)
  rownames(out) <- NULL
  out <- out[order(out[[id]], out[[start]], out[[stop]], out$.tcoarsen_orig_row),
             , drop = FALSE]

  if (length(absorb_vars) > 0L) {
    for (v in absorb_vars) {
      out[[v]] <- stats::ave(as.integer(out[[v]]), out[[id]], FUN = function(x) cummax(x))
    }
  }

  if (any(out[[stop]] <= out[[start]])) {
    stop("Coarsening produced rows with `stop <= start`. This should not happen; ",
         "please inspect the input data and coarsening settings.", call. = FALSE)
  }

  if (isTRUE(add_visit)) {
    out[[visit_name]] <- as.integer(floor((out[[start]] - origin) / bin_width +
                                            sqrt(.Machine$double.eps)))
  }

  out$.tcoarsen_orig_row <- NULL
  out$.tcoarsen_snap_start <- NULL

  tcoarsen_info <- NULL
  if (isTRUE(diagnostics)) {
    tcoarsen_info <- list(
      call = match.call(),
      id = id,
      start = start,
      stop = stop,
      event = event,
      covs = covs,
      bin_width = bin_width,
      dir_coarsen = dir_coarsen,
      origin = origin,
      absorb_vars = absorb_vars,
      keep_terminal_time = keep_terminal_time,
      n_subjects = length(unique(out[[id]])),
      n_rows_original = n_original,
      n_rows_returned = nrow(out),
      n_records_after_bin_collapse = diag_env$n_records_after_bin_collapse,
      n_inserted_lkcf_rows = max(0L, nrow(out) - diag_env$n_records_after_bin_collapse),
      n_collapsed_rows = diag_env$n_collapsed_rows,
      n_collapsed_bins = diag_env$n_collapsed_bins,
      n_dropped_updates_after_exit = diag_env$n_dropped_updates_after_exit,
      n_subjects_with_post_exit_updates = diag_env$n_subjects_with_post_exit_updates,
      visit_summary = if (isTRUE(add_visit)) {
        .tcoarsen_visit_summary(out, id = id, event = event,
                                visit_name = visit_name,
                                absorb_vars = absorb_vars)
      } else {
        NULL
      }
    )
    #attr(out, "tcoarsen") <- tcoarsen_info
  }

  if (isTRUE(verbose)) {
    msg <- paste0(
      "tcoarsen(): returned ", nrow(out), " rows for ",
      length(unique(out[[id]])), " subjects; bin_width = ",
      format(bin_width, trim = TRUE, scientific = FALSE),
      ", dir_coarsen = '", dir_coarsen, "'."
    )
    message(msg)
    if (diag_env$n_collapsed_rows > 0L) {
      message("tcoarsen(): ", diag_env$n_collapsed_rows,
              " row(s) collapsed because multiple updates mapped to the same grid time.")
    }
    if (diag_env$n_dropped_updates_after_exit > 0L) {
      message("tcoarsen(): ", diag_env$n_dropped_updates_after_exit,
              " coarsened update row(s) were at or after the terminal exit time and were dropped.")
    }
  }
  out <- list(dat_coarsen = out, diagnostics = tcoarsen_info)
  class(out) <- "tcoarsen"
  out
}

.tcoarsen_check_string <- function(x, nm) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
    stop("`", nm, "` must be a single non-missing character string.",
         call. = FALSE)
  }
  invisible(TRUE)
}

.tcoarsen_to_numeric_time <- function(x, nm) {
  if (inherits(x, "Date") || inherits(x, "POSIXt")) {
    stop("`", nm, "` must be numeric. Convert Date/POSIX times to a numeric ",
         "time scale before calling `tcoarsen()`.", call. = FALSE)
  }
  z <- suppressWarnings(as.numeric(x))
  if (anyNA(z)) {
    stop("`", nm, "` could not be converted to a non-missing numeric vector.",
         call. = FALSE)
  }
  z
}

.tcoarsen_to_binary01 <- function(x, nm, allow_na = FALSE) {
  if (is.logical(x)) {
    z <- as.integer(x)
  } else if (is.factor(x)) {
    z <- suppressWarnings(as.integer(as.character(x)))
  } else if (is.character(x)) {
    z <- suppressWarnings(as.integer(trimws(x)))
  } else if (is.numeric(x) || is.integer(x)) {
    z <- as.integer(x)
  } else {
    stop("`", nm, "` must be coded as 0/1, logical, or a factor/character ",
         "with levels/values '0' and '1'.", call. = FALSE)
  }

  if (!allow_na && anyNA(z)) {
    stop("`", nm, "` cannot contain missing values and must be coded 0/1.",
         call. = FALSE)
  }
  ok <- is.na(z) | z %in% c(0L, 1L)
  if (!all(ok)) {
    stop("`", nm, "` must be coded as 0/1.", call. = FALSE)
  }
  z
}

.tcoarsen_check_intervals <- function(df,
                                      id,
                                      start,
                                      stop,
                                      gap_action,
                                      tol) {
  split_ids <- split(df, df[[id]], drop = TRUE)
  ids_with_gaps <- character(0L)
  ids_with_overlaps <- character(0L)

  for (nm in names(split_ids)) {
    d <- split_ids[[nm]]
    d <- d[order(d[[start]], d[[stop]], d$.tcoarsen_orig_row), , drop = FALSE]
    if (nrow(d) <= 1L) next
    previous_stop <- d[[stop]][seq_len(nrow(d) - 1L)]
    next_start <- d[[start]][seq.int(2L, nrow(d))]
    if (any(next_start < previous_stop - tol)) {
      ids_with_overlaps <- c(ids_with_overlaps, nm)
    }
    if (any(next_start > previous_stop + tol)) {
      ids_with_gaps <- c(ids_with_gaps, nm)
    }
  }

  if (length(ids_with_overlaps) > 0L) {
    stop("Found overlapping intervals for subject(s): ",
         paste(utils::head(ids_with_overlaps, 5L), collapse = ", "),
         if (length(ids_with_overlaps) > 5L) " ..." else "",
         ". Overlapping start-stop intervals are ambiguous for LKCF coarsening.",
         call. = FALSE)
  }

  if (length(ids_with_gaps) > 0L && gap_action != "ignore") {
    msg <- paste0(
      "Found gaps between adjacent start-stop intervals for subject(s): ",
      paste(utils::head(ids_with_gaps, 5L), collapse = ", "),
      if (length(ids_with_gaps) > 5L) " ..." else "",
      ". `tcoarsen()` assumes continuous at-risk follow-up when carrying values forward."
    )
    if (gap_action == "stop") {
      stop(msg, call. = FALSE)
    } else if (gap_action == "warn") {
      warning(msg, call. = FALSE)
    }
  }

  invisible(TRUE)
}

.tcoarsen_visit_summary <- function(out,
                                    id,
                                    event,
                                    visit_name,
                                    absorb_vars = NULL) {
  by_visit <- data.frame(
    visit = out[[visit_name]],
    rows = rep(1L, nrow(out)),
    ids = out[[id]],
    events = out[[event]]
  )

  row_event <- stats::aggregate(cbind(rows, events) ~ visit,
                                data = by_visit,
                                FUN = sum)
  id_count <- stats::aggregate(ids ~ visit,
                               data = by_visit,
                               FUN = function(z) length(unique(z)))
  names(id_count)[names(id_count) == "ids"] <- "n_ids"

  ans <- merge(row_event, id_count, by = "visit", all = TRUE)
  ans <- ans[order(ans$visit), , drop = FALSE]

  if (length(absorb_vars) > 0L) {
    for (v in absorb_vars) {
      if (!v %in% names(out)) next
      active_df <- data.frame(
        visit = out[[visit_name]],
        active = as.integer(out[[v]] == 1L)
      )
      active_sum <- stats::aggregate(active ~ visit,
                                     data = active_df,
                                     FUN = sum)
      names(active_sum)[names(active_sum) == "active"] <- paste0(v, "_active_rows")
      ans <- merge(ans, active_sum, by = "visit", all = TRUE)

      incident <- stats::ave(as.integer(out[[v]]), out[[id]], FUN = function(z) {
        z <- as.integer(z)
        if (length(z) == 0L) return(integer(0L))
        as.integer(c(z[1L] == 1L, diff(z) == 1L))
      })
      incident_df <- data.frame(
        visit = out[[visit_name]],
        incident = as.integer(incident)
      )
      incident_sum <- stats::aggregate(incident ~ visit,
                                       data = incident_df,
                                       FUN = sum)
      names(incident_sum)[names(incident_sum) == "incident"] <- paste0(v, "_incident_rows")
      ans <- merge(ans, incident_sum, by = "visit", all = TRUE)
    }
    ans <- ans[order(ans$visit), , drop = FALSE]
  }

  rownames(ans) <- NULL
  ans
}
