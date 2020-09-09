#' Suppress function messages and Concatenate and Print (cat)
#'
#' This function is used to suppress information printed from external functions
#' that make internal use of \code{link{message}} and \code{\link{cat}}, which
#' provide information in interactive R sessions. For simulations, the session
#' is not interactive, and therefore this type of output should be suppressed.
#' For similar behavior for suppressing warning messages see
#' \code{\link{suppressWarnings}}, though use this function carefully as some
#' warnings can be meaningful and unexpected.
#' @name quiet
#' @aliases quiet
#' @param ... the functional expression to be evaluated
#' @param messages logical; suppress all messages?
#' @param cat logical; suppress all concatenate and print calls from \code{\link{cat}}?
#' @author Phil Chalmers
#' @return quiet
#' @keywords internal
#' @references Sigal, M. J., & Chalmers, R. P. (2016). Play it again: Teaching statistics with Monte Carlo simulation. \code{Journal of Statistics Education, 24}(3), 136-156.
quiet <- function(..., messages=FALSE, cat=FALSE){
  if(!cat){
    sink(tempfile())
    on.exit(sink())
  }
  out <- if(messages) eval(...) else suppressMessages(eval(...))
  out
}
