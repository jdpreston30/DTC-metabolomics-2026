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

#' Round to 1 Decimal Place (Half Up)
#'
#' Rounds numbers to 1 decimal place using the "round half up" rule.
#' Values where the second decimal place is 5 or greater round up,
#' values less than 5 round down.
#'
#' @param x Numeric vector to round
#'
#' @return Numeric vector with values rounded to 1 decimal place
#'
#' @details
#' This function ensures proper rounding to 1 decimal place:
#' - 1.45 rounds to 1.5
#' - 1.44 rounds to 1.4
#' - 5.55 rounds to 5.6
#' - 2.84 rounds to 2.8
#'
#' Handles floating point precision issues by first rounding to 2 decimal places
#' to clean up floating point errors, then applies half-up rounding to 1 decimal.
#'
#' @examples
#' round_1dp(c(1.45, 1.44, 5.55, 2.84)) # Returns: 1.5 1.4 5.6 2.8
#'
#' @export
round_1dp <- function(x) {
  # First round to 2 decimal places to handle floating point errors
  x_2dp <- round(x, 2)
  # Then apply half-up rounding to 1 decimal place
  round_half_up(x_2dp * 10) / 10
}
