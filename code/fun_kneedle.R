#' Detect Knee Point in a Curve Using the Kneedle Algorithm
#'
#' This function implements the Kneedle algorithm to detect the "knee" or "elbow"
#' point in a curve, as described in Satopaa et al. (2011). The knee point is the
#' point of maximum curvature, often used to identify a transition in data behavior
#' (e.g., optimal clustering parameters). The algorithm normalizes the input data,
#' computes the difference between the curve and a reference line, and identifies
#' the knee based on peaks in the difference curve, adjusted by a sensitivity parameter.
#'
#' The code is adapted from \url{https://github.com/etam4260/kneedle}.
#'
#' @param x A numeric vector of x coordinates, strictly increasing or decreasing.
#' @param y A numeric vector of y coordinates, same length as `x`.
#' @param decreasing Logical, indicating if the curve is decreasing (TRUE) or
#'   increasing (FALSE). If NULL, the function estimates the direction based on the
#'   difference between the first and last y values (default: NULL).
#' @param concave Logical, indicating if the curve is concave (TRUE) or convex
#'   (FALSE). If NULL, the function estimates concavity based on the average second
#'   derivative of the data (default: NULL).
#' @param sensitivity Numeric, sensitivity for detecting the knee point. Higher
#'   values make the algorithm more selective, requiring a stronger peak to identify
#'   the knee (default: 1).
#' @return A numeric vector of length 2 containing the x and y coordinates of the
#'   detected knee point. If no knee is found, returns `c(NA, NA)`.
#' @details
#' The Kneedle algorithm, described in Satopaa et al. (2011), detects the knee point
#' by normalizing the input data to [0, 1], computing the difference between the
#' normalized curve and a reference line (x = y or x = 1 - y, depending on
#' direction and concavity), and identifying peaks in the difference curve. The first
#' peak exceeding a threshold (adjusted by `sensitivity`) is selected as the knee.
#' The function automatically estimates `decreasing` and `concave` if not specified,
#' using the slope from the first to last point and the average second derivative,
#' respectively.
#'
#' @references
#' Satopaa, V., Albrecht, J., Irwin, D., & Raghavan, B. (2011). Finding a "Kneedle"
#' in a Haystack: Detecting Knee Points in System Behavior. *31st International
#' Conference on Distributed Computing Systems Workshops*, 166-171.
#' \doi{10.1109/ICDCSW.2011.20}
#'
#' @examples
#' \dontrun{
#' # Example with synthetic data
#' x <- 1:10
#' y <- c(1, 2, 3, 4, 10, 20, 30, 40, 50, 60)
#' knee <- kneedle(x, y, sensitivity = 1)
#' plot(x, y, type = "l", main = "Knee Point Detection")
#' points(knee[1], knee[2], col = "red", pch = 19)
#' }
#' @importFrom quantmod findPeaks
#' @importFrom stats diff
#' @export
#'
kneedle <- function(x, y, decreasing = NULL, concave = NULL, sensitivity = 1) {
    # Input validation
    if (!is.numeric(x) || !is.numeric(y)) {
        stop("Both `x` and `y` must be numeric vectors.")
    }
    if (length(x) == 0 || length(y) == 0) {
        stop("Both `x` and `y` must have length greater than 0.")
    }
    if (length(x) != length(y)) {
        stop("`x` and `y` must have the same length.")
    }
    if (any(!is.finite(x) | !is.finite(y))) {
        stop("`x` and `y` must not contain NA, NaN, or Inf values.")
    }
    if (!is.numeric(sensitivity) || sensitivity <= 0) {
        stop("`sensitivity` must be a positive numeric value.")
    }
    if (!requireNamespace("quantmod", quietly = TRUE)) {
        stop("Package 'quantmod' is required. Please install it using install.packages('quantmod').")
    }

    # Combine and sort data by x
    data <- matrix(c(x, y), ncol = 2)
    data <- data[order(data[, 1], decreasing = FALSE), ]

    # Estimate decreasing if not specified
    if (is.null(decreasing)) {
        decreasing <- data[nrow(data), 2] < data[1, 2]
        # message("`decreasing` not specified. Estimated as ", decreasing,
        #         " based on first and last y values.")
    }

    # Estimate concave if not specified
    if (is.null(concave)) {
        second_deriv <- diff(diff(data[, 2]) / diff(data[, 1]))
        concave <- mean(second_deriv, na.rm = TRUE) > 0
        # message("`concave` not specified. Estimated as ", concave,
        #         " based on average second derivative.")
    }

    # Normalize x and y to [0,1]
    max_x <- max(data[, 1])
    min_x <- min(data[, 1])
    max_y <- max(data[, 2])
    min_y <- min(data[, 2])
    if (max_x == min_x || max_y == min_y) {
        warning("`x` or `y` has no variation. Returning NA for knee point.")
        return(c(NA, NA))
    }
    data[, 1] <- (data[, 1] - min_x) / (max_x - min_x)
    data[, 2] <- (data[, 2] - min_y) / (max_y - min_y)

    # Compute difference curve based on direction and concavity
    if (concave && !decreasing) {
        differ <- abs(data[, 2] - data[, 1])
    } else if (concave && decreasing) {
        differ <- abs(data[, 2] - (1 - data[, 1]))
    } else if (!concave && !decreasing) {
        differ <- abs(data[, 2] - data[, 1])
    } else if (!concave && decreasing) {
        differ <- abs(data[, 2] - (1 - data[, 1]))
    }
    data <- cbind(data, differ)

    # Find peaks in the difference curve
    peak_indices <- quantmod::findPeaks(differ) - 1
    if (length(peak_indices) == 0) {
        warning("No peaks found in difference curve. Returning NA for knee point.")
        return(c(NA, NA))
    }

    # Compute threshold for knee detection
    diff_x <- diff(data[, 1])
    T_lm_x_s <- sensitivity * mean(diff_x, na.rm = TRUE)

    # Identify knee point
    knee <- NULL
    for (i in seq_along(peak_indices)) {
        Ti <- data[peak_indices[i], 3] - T_lm_x_s
        for (j in (peak_indices[i]):if (i + 1 < length(peak_indices)) peak_indices[i + 1] else length(differ)) {
            if (differ[j] < Ti) {
                knee <- peak_indices[i]
                break
            }
        }
        if (!is.null(knee)) break
    }

    # Return NA if no knee is found
    if (is.null(knee)) {
        warning("No knee point identified. Returning NA.")
        return(c(NA, NA))
    }

    # Denormalize and return knee point coordinates
    x_knee <- (max_x - min_x) * data[knee, 1] + min_x
    y_knee <- (max_y - min_y) * data[knee, 2] + min_y
    return(c(x = as.numeric(x_knee), y = as.numeric(y_knee)))
}
