#' Round Half Up (Banker's Rounding Alternative)
#'
#' Rounds numbers using the "round half up" rule where values exactly halfway
#' between two integers are always rounded up (toward positive infinity).
#' This is an alternative to R's default "round half to even" behavior.
#'
#' @param x Numeric vector to round
#'
#' @return Numeric vector with values rounded using the "round half up" rule
#'
#' @details
#' R's default rounding uses "round half to even" (banker's rounding), where
#' 0.5 rounds to 0 and 1.5 rounds to 2. This function implements "round half up"
#' where 0.5 rounds to 1 and 1.5 rounds to 2.
#'
#' @examples
#' # Compare with R's default rounding
#' round(c(0.5, 1.5, 2.5)) # Returns: 0 2 2 (round half to even)
#' round_half_up(c(0.5, 1.5, 2.5)) # Returns: 1 2 3 (round half up)
#'
#' @export
round_half_up <- function(x) {
  floor(x + 0.5)
}
