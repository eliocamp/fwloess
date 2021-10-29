#' Local regression with fixed width
#'
#' [lowess] it's awesome, but the span parameter is hard to interpret as the
#' actual points that are used depends on the size of the series. Furthermore,
#' the actual distance each point is to the centre of the regression is not
#' taken into account, which is a problem for series irregularly sampled.
#'
#' @param formula a formula specifying the numeric response and one numeric
#' predictor.
#' @param data an optional data frame.
#' @param weights currently ignored.
#' @param weighting whether to use bicubic weighting or not.
#' @param span width of the smoothing window.
#' @param degree degree of the polynomial.
#' @param min_n minimum number of points needed to perform local regression.
#'
#'
#'
#' @export
fw_loess <- function(data = NULL, formula, span, degree = 1, min_n = degree + 3, weights, weighting = TRUE) {

  data <- stats::model.frame(formula, data = data)
  data <- data[stats::complete.cases(data), ]

  vars <- all.vars(formula)
  y_var <- vars[1]
  x_var <- vars[2]

  span_seq <- seq(min(data[[x_var]]),
                  max(data[[x_var]]),
                  by = span)
  span <- as.numeric(diff(span_seq)[1])

  structure(.Data = list(data = data,
                         span = span,
                         degree = degree,
                         min_n = min_n,
                         weighting = weighting,
                         vars = vars,
                         y_var = y_var,
                         x_var = x_var), class = c("fw_loess"))

}

