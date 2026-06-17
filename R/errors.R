# SPDX-FileCopyrightText: The fastMethyl authors
# SPDX-License-Identifier: Artistic-2.0
#
# Centralised error helpers so every failure a user can hit surfaces with the
# same shape: what the package was doing, the underlying cause, what to try, and
# -- when it looks like a genuine bug -- how to report it.
#
# Two condition classes carry the distinction analyze()'s phase wrappers need:
#
#   fastMethyl_user_error -- an actionable problem with the caller's inputs,
#     files, permissions, environment, or arguments. The message already says
#     exactly what to fix, so the phase wrappers relay it verbatim and never bury
#     it under a "this looks like a bug" footer.
#
#   fastMethyl_error -- an error caught by a phase wrapper (.withPhase) and
#     augmented with the phase context, an optional remediation hint, and the
#     bug-report footer. This is what an *unexpected* failure (anything thrown by
#     minfi / gdsfmt / bigmelon / base R that we did not raise ourselves as a
#     user error) is turned into.

# Where to report a suspected bug. Kept in sync with DESCRIPTION's BugReports.
.bugReportUrl <- "https://github.com/sebhoss/fastMethyl/issues"

# Raise an actionable, user-facing error. `...` are pasted into the message
# (like stop()). Tagged fastMethyl_user_error so .withPhase passes it through
# unchanged -- its message is already the remedy, so it must not be wrapped.
.userStop <- function(...) {
  stop(structure(
    class = c("fastMethyl_user_error", "error", "condition"),
    list(message = paste0(...), call = NULL)
  ))
}

# sprintf-formatting convenience: .userStopf(fmt, ...) is
# .userStop(sprintf(fmt, ...)). Keeps value-bearing messages on one (flat) line.
.userStopf <- function(fmt, ...) {
  .userStop(sprintf(fmt, ...))
}

# Assert `cond`; on failure raise `msg` as a user error. The classed equivalent
# of a single stopifnot() expression, so a violated precondition reads as the
# actionable problem it is rather than as an internal assertion.
.check <- function(cond, ...) {
  if (!isTRUE(cond)) {
    .userStop(...)
  }
  invisible(TRUE)
}

# Build (but do not raise) the augmented message for an unexpected failure in a
# given phase. Separated out so it can be unit-tested without a running error.
.phaseMessage <- function(doing, cause, hint = NULL) {
  parts <- c(
    sprintf("fastMethyl could not complete while %s.", doing),
    "",
    "  What went wrong:",
    paste0("    ", cause)
  )
  if (!is.null(hint) && nzchar(hint)) {
    parts <- c(parts, "", "  What you can try:", paste0("    ", hint))
  }
  parts <- c(
    parts, "",
    paste0(
      "If this looks like a bug in fastMethyl rather than a problem with your ",
      "inputs, please report it at ", .bugReportUrl, " -- include the error ",
      "above and the output of utils::sessionInfo()."
    )
  )
  paste(parts, collapse = "\n")
}

# Run `expr`; if it raises an *unexpected* error, abort with a message naming
# the phase (`doing`, an "-ing" clause), the underlying cause, and `hint`. A
# fastMethyl_user_error (already actionable) and a fastMethyl_error (already
# wrapped by an inner phase) are re-raised unchanged so context is added exactly
# once and user errors are never relabelled as bugs.
#
# The error is captured into a value and re-raised from the function body, NOT
# from inside the tryCatch handler: re-signalling with stop() *inside* a handler
# is caught by that same tryCatch's other handlers, which would re-wrap a user
# error as an unexpected one. Returning to the body first sidesteps that.
.withPhase <- function(doing, expr, hint = NULL) {
  caught <- tryCatch(
    list(ok = TRUE, value = expr),
    error = function(e) list(ok = FALSE, error = e)
  )
  if (caught$ok) {
    return(caught$value)
  }
  e <- caught$error
  # Our own conditions are already actionable / already wrapped: pass through.
  if (inherits(e, "fastMethyl_user_error") || inherits(e, "fastMethyl_error")) {
    stop(e)
  }
  stop(structure(
    class = c("fastMethyl_error", "error", "condition"),
    list(message = .phaseMessage(doing, conditionMessage(e), hint), call = NULL)
  ))
}
