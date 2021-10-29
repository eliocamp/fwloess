
#' @export
predict.fw_loess <- function(fit, newdata = NULL, se = FALSE, ...) {
  if (is.null(newdata)) {
    x.out_org <- fit$data[[fit$x_var]]
    x.out <- as.numeric(x.out_org)
  } else if (is.character(newdata)) {
    x.out_org <- seq(min(fit$data[[fit$x_var]]),
                     max(fit$data[[fit$x_var]]),
                     by = newdata)
    x.out <- as.numeric(x.out_org)
  } else {
    x.out_org <- newdata
    x.out <- as.numeric(newdata)
  }

  fit$data[[fit$x_var]] <- as.numeric(fit$data[[fit$x_var]])
  df <- sd.out <- rep(NA, length = length(x.out))
  y.out <- matrix(NA, ncol = fit$degree + 1, nrow = length(x.out))
  half_span <- fit$span/2
  # min_n <- degree + 2

  for (i in seq_along(x.out)) {
    sub <- fit$data[[fit$x_var]] >= (x.out[i] - half_span) & fit$data[[fit$x_var]] <= (x.out[i] + half_span)
    if (sum(sub) < fit$min_n) {
      y.out[i] <- NA
    } else{
      # ss <<- sub

      if (isTRUE(fit$weighting)) {
        weight_x <- (1 - abs(fit$data[[fit$x_var]][sub] - x.out[i])^3/half_span^3)^3
        weight <- weight_x
      } else {
        weight <- rep(1, sum(sub))
      }


      x_min <- min(fit$data[[fit$x_var]][sub])

      # X <- poly(x.out - x_min, degree = degree, raw = TRUE)
      X <- vapply(c(0, seq_len(fit$degree)), function(d) (fit$data[[fit$x_var]][sub] - x_min)^d,
                  numeric(length = length(fit$data[[fit$x_var]][sub])))

      # X.out <- poly(x.out[i] - x_min, degree = degree, raw = TRUE)
      X.out <- vapply(c(0, seq_len(fit$degree)), function(d) (x.out[i] - x_min)^d,
                      numeric(length = 1))

      local_fit <- stats::lm.wfit(x = X,
                                  y =  fit$data[[fit$y_var]][sub],
                                  w = weight)

      y.out[i, 1] <- local_fit$coefficients %*% X.out


      # derivates
      y.out[i, -1] <- vapply(seq_len(fit$degree), function(d) {
        I <- seq_along(local_fit$coefficients) - 1
        a <- local_fit$coefficients[I >= d]
        I <- I[I >= d]


        b <- factorial(I)/factorial(I - d)
        a %*% (b* (x.out[i] - x_min)^(I - d))
      }, 1)

      df[i] <- sum(sub) - fit$degree - 1
      sd.out[i] <- stats::sd(local_fit$residuals)*sqrt(t(X.out) %*% solve(t(X) %*% X) %*% X.out)
    }
  }
  # browser()

  out <- as.data.frame(y.out)
  if (fit$degree != 0) {
    names_derv <- paste0(fit$vars[1], "_d", seq_len(fit$degree), fit$vars[2])
  } else {
    names_derv <- NULL
  }

  colnames(out) <- c(fit$vars[1], names_derv)
  x <- data.frame(x = x.out_org)
  colnames(x) <- fit$vars[2]
  errors <- data.frame(se = sd.out)

  out <- cbind(x, out, errors, df = df)

  return(out)
}


#' @export
predictdf.fw_loess <- function(model, xseq, se, level) {
  pred <- stats::predict(model, newdata = xseq)

  if (se == TRUE) {
    ci <- pred$se * stats::qt(level/2 + 0.5, pred$df)
    ymin = pred$y - ci
    ymax = pred$y + ci
    pred$ymin <- ymin
    pred$ymax <- ymax

  }
  pred

}


#' @export
fitted.fw_loess <- function(object, ...) {
  stats::predict(object)[[object$y_var]]
}
